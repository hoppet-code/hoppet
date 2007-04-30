module new_as
  use types; use consts_dp
  use assertions; use warnings_and_errors
  implicit none
  private

  type na_segment
     real(dp) :: tlo, thi, dt
     real(dp) :: beta0, beta1, beta2
     real(dp), pointer :: ra(:)
     !-- above iflip we evolve upwards, below iflip we evolve downwards
     integer  :: iflip, dummy
  end type na_segment
  
  type na_handle
     private
     type(na_segment), pointer :: seg(:)
     integer                   :: nlo, nhi
     integer                   :: nloop, fixnf
     !--- these are startoff parameters
     real(dp)                  :: alfas, Q
     real(dp)                  :: quark_masses(1:6)
     real(dp)                  :: muMatch_mQuark
  end type na_handle
  
  !-- this seems to give reasonable results (c. 10^-9 precision
  !   almost everywhere if we start at Z peak)
  real(dp) :: dt_base = 0.2_dp
  !real(dp), parameter :: dt_base = 0.1_dp
  !--------- go from 0.5 GeV to over 10^20 GeV
  real(dp), parameter :: tlo = -1.3862944 !-3.2188758_dp
  real(dp), parameter :: thi = 93.0_dp 
  integer,  parameter :: nofixnf = -1000000045

  public :: na_handle
  public :: na_Init, na_Value, na_Del, na_NumberOfLoops
  public :: na_nfRange, na_nfAtQ, na_QrangeAtNf, na_QuarkMass
  public :: na_Set_dt_base

  type(na_segment), pointer :: seg

  interface SetDefaultCouplingDt
     module procedure na_Set_dt_base
  end interface
  public :: SetDefaultCouplingDt, DefaultCouplingDt
     

contains
  
  !======================================================================
  !! Set the value of dt_base (used to decide spacings for new coupling
  !! instances).
  subroutine na_Set_dt_base(new_dt_base)
    real(dp) :: new_dt_base
    dt_base = new_dt_base
  end subroutine na_Set_dt_base
  


  !---------------------------------------------------
  !! Set the value of dt_base (used to decide spacings for new coupling
  !! instances).
  function DefaultCouplingDt() 
    real(dp) :: DefaultCouplingDt
    DefaultCouplingDt = dt_base
  end function DefaultCouplingDt
  

  !---------------------------------------------------
  subroutine na_Init(nah, alfas, Q, nloop, fixnf, quark_masses, muMatch_mQuark)
    use qcd; 
    type(na_handle),  intent(out), target :: nah
    real(dp), intent(in), optional :: alfas, Q
    integer,  intent(in), optional :: nloop, fixnf
    real(dp), intent(in), optional :: quark_masses(4:6)
    real(dp), intent(in), optional :: muMatch_mQuark
    !---------------------------------------------
    !real(dp) :: alfas_lcl, Q_lcl
    !integer  :: nloop_lcl, fixnf_lcl
    !-------------------------------
    real(dp) :: alfas_here, lnmatch
    real(dp) :: tstart, t, ra
    integer  :: nf_store
    integer  :: nbin, i, j, nseg, istart
    integer  :: nbin_total

    !-- we may well play with nf, so need to be able to reset it
    nf_store = nf_int

    nah%alfas = default_or_opt(0.118_dp, alfas)
    nah%Q     = default_or_opt(91.2_dp, Q)
    nah%nloop = default_or_opt(2, nloop)
    
    nah%fixnf = default_or_opt(nofixnf, fixnf)

    nah%muMatch_mQuark = default_or_opt(one, muMatch_mQuark)

    nah%quark_masses(1:3) = quark_masses_def(1:3)
    if (present(quark_masses)) then
       nah%quark_masses(4:6) = quark_masses(4:6)
    else
       nah%quark_masses(4:6) = quark_masses_def(4:6)
    end if
    

    !-- work out which ranges are needed
    !   even in the fixnf case these are used as guides
    !   for when to switch dt?
    do i = lbound(nah%quark_masses, dim=1), ubound(nah%quark_masses, dim=1)
       if (tlo > tOfQ(nah%muMatch_mQuark*nah%quark_masses(i))) nah%nlo = i
       if (thi > tOfQ(nah%muMatch_mQuark*nah%quark_masses(i))) nah%nhi = i
    end do
    
    !-- now set up the parameters of the segments
    allocate(nah%seg(nah%nlo:nah%nhi))
    nbin_total = 0
    do i = nah%nlo, nah%nhi
       seg => nah%seg(i)
       nah%seg(i)%tlo = max(tlo, tOfQ(nah%muMatch_mQuark*nah%quark_masses(i)))
       if (i < nah%nhi) then
          seg%thi = tOfQ(nah%muMatch_mQuark*nah%quark_masses(i+1))
       else
          seg%thi = thi
       end if

       !-- now fix running coupling
       if (nah%fixnf == nofixnf) then
          call qcd_SetNf(i)
       else
          call qcd_SetNf(nah%fixnf)
       end if
       
       !-- lazy person's option for getting number of loops automatically
       !   in derivative calculation
       seg%beta0 = beta0
       select case(nah%nloop)
       case(1)
          seg%beta1 = zero
          seg%beta2 = zero
       case(2)
          seg%beta1 = beta1
          seg%beta2 = zero
       case(3)
          seg%beta1 = beta1
          seg%beta2 = beta2
       case default
          call wae_Error('na_init: unsupported number of loops requested')
       end select
       
       !-- first guess, which then needs to be "adapted" to the range
       seg%dt    = dt_base * half**(6-i)
       nbin = ceiling((seg%thi - seg%tlo)/seg%dt)
       nbin_total = nbin_total + nbin
       seg%dt    = (seg%thi - seg%tlo) / nbin
       allocate(seg%ra(0:nbin))
    end do
    !write(0,*) 'new_as: used total bins ', nbin_total

    !-- find out in which segment we start
    tstart = tOfQ(nah%Q)
    do i = nah%nlo, nah%nhi
       if (tstart <= nah%seg(i)%thi .and. tstart >= nah%seg(i)%tlo) exit
    end do
    if (i > nah%nhi) &
         &call wae_Error('na_init: Specified Q in not in supported range')
    nseg = i

    !-- fill up the starting segment
    ra = one/nah%alfas
    seg => nah%seg(nseg)
    istart = nint((tstart-seg%tlo)/seg%dt)
    t = tstart
    seg%ra(istart) = na_evolve(ra, seg%tlo+istart*seg%dt - tstart)
    do i = istart+1, ubound(seg%ra,dim=1)
       seg%ra(i) = na_evolve(seg%ra(i-1), seg%dt)
    end do
    do i = istart-1, lbound(seg%ra,dim=1), -1
       seg%ra(i) = na_evolve(seg%ra(i+1), -seg%dt)
    end do
    seg%iflip = istart

    !-- fill up the segments above
    do j = nseg+1, nah%nhi
       seg => nah%seg(j)
       alfas_here = one/nah%seg(j-1)%ra(ubound(nah%seg(j-1)%ra,dim=1))
       if (mass_steps_on .and. nah%fixnf == nofixnf) then
          !write(0,*) '------ changing as from ',j-1,' to ',j
          !write(0,*) j-1, alfas_here
          lnmatch = two*log(nah%muMatch_mQuark)
          select case (nah%nloop)
          case(2)
             alfas_here = alfas_here*(1 + alphastep11 * lnmatch * alfas_here)
          case(3)
             alfas_here = alfas_here*(1 + alphastep11 * lnmatch * alfas_here&
                  & + (alphastep22 * lnmatch**2 + alphastep21 * lnmatch&
                  &    + alphastep20_pole)*alfas_here**2)
             !alfas_here = alfas_here*(1 + alphastep20_pole*alfas_here**2)
          end select
          !write(0,*) j, alfas_here
       end if
       seg%ra(0) = one/alfas_here
       ! recall that this is the reciprocal of alpha!
       do i = 1, ubound(seg%ra,dim=1)
          seg%ra(i) = na_evolve(seg%ra(i-1), seg%dt)
       end do
       seg%iflip = -1
    end do
    
    !-- fill up the segments below
    do j = nseg-1, nah%nlo, -1
       seg => nah%seg(j)
       alfas_here = one/nah%seg(j+1)%ra(0)
       if (mass_steps_on .and. nah%fixnf == nofixnf) then
          !write(0,*) '------ changing as from ',j+1,' to ',j
          !write(0,*) j+1, alfas_here
          lnmatch = two*log(nah%muMatch_mQuark)
          select case (nah%nloop)
          case(2)
             alfas_here = alfas_here*(1 - alphastep11 * lnmatch * alfas_here)
          case(3)
             alfas_here = alfas_here*(1 - alphastep11 * lnmatch * alfas_here&
                  & + (alphastep22 * lnmatch**2 - alphastep21 * lnmatch&
                  &    - alphastep20_pole)*alfas_here**2)
             !alfas_here = alfas_here*(1 - alphastep20_pole*alfas_here**2)
          end select
          !write(0,*) j, alfas_here
       end if
       seg%ra(ubound(seg%ra,dim=1)) = one/alfas_here
       do i = ubound(seg%ra,dim=1)-1, 0, -1
          seg%ra(i) = na_evolve(seg%ra(i+1), -seg%dt)
       end do
       seg%iflip = ubound(seg%ra,dim=1)+1
    end do

    !--- set things back to their origins
    call qcd_SetNf(nf_store)
  end subroutine na_init
  

  !-------------------------------------------------------------------------
  ! If fixnf is present, then the program will force (within limits of +-dt)
  ! the alpha that is returned to be that corresponding to fixnf flavours.
  ! 
  ! With smooth matching conditions, this is not a big deal, but with
  ! steps at alpha_s thresholds, it is important so that certain D.E.
  ! routines (e.g. for PDF evolution) do not give alpha values at end 
  ! points that give non-smoothness.
  function na_Value(nah, Q, fixnf) result(res)
    type(na_handle), intent(in), target :: nah
    real(dp),        intent(in)         :: Q
    integer,         intent(in), optional :: fixnf
    real(dp)                            :: res
    !-------------------------------------
    real(dp) :: t, trescaled, ra, delta_t
    integer  :: nseg, i, n
    integer, save :: warn_id = warn_id_INIT
    integer, parameter :: max_warn = 1

    t = tOfQ(Q)

    if (present(fixnf) .and. nah%fixnf == nofixnf) then
       if (fixnf < nah%nlo .or. fixnf > nah%nhi) then
          call wae_error('na_Value:', 'the fixnf requested is&
               & outside the range supported this na_handle')
       end if
       if (t < tlo .or. t > thi) then
          call wae_error('na_Value:', 'the Q value is&
               & outside the range supported this na_handle')
       end if
       nseg = fixnf
       !if (t > nah%seg(nseg)%thi+ nah%seg(nseg)%dt .or.&
       !     & t < nah%seg(nseg)%tlo-nah%seg(nseg)%dt) then
       !   call wae_error('na_Value:', &
       !        &  'With fixnf, Q is too far outside supported range.')
       !end if
    else
       if (present(fixnf) .and. nah%fixnf /= nofixnf) then
          if (fixnf /= nah%fixnf) then
             call wae_error('na_Value:', 'the fixnf requested is &
                  &different from that supported by na_handle')
          end if
       end if
       do nseg = nah%nlo, nah%nhi
          if (t <= nah%seg(nseg)%thi .and. t >= nah%seg(nseg)%tlo) exit
       end do
       if (nseg > nah%nhi) &
            &call wae_Error('na_Value: Specified Q is not in supported range'&
            &,dbleval=Q)
    end if
    
    seg => nah%seg(nseg)
    trescaled = (t - seg%tlo)/seg%dt
    if (trescaled > seg%iflip) then
       i = floor(trescaled)
    else
       i = ceiling(trescaled)
    end if
    i = max(0,min(ubound(seg%ra,dim=1),i))

    !-- support going well beyond supported limits of this nf, even
    !   if the procedure is not particularly recommended.
    delta_t = t - (seg%tlo+i*seg%dt)
    if (abs(delta_t) <= 1.3_dp*seg%dt) then
       res = one/na_evolve(seg%ra(i), delta_t)
    else
       call wae_warn(max_warn,warn_id,'na_Value: will evolve &
            &fixed-nf alpha_s beyond precalculated range.',&
            &'This procedure may be very slow')
       !write(0,*) Qoft(seg%tlo),Qoft(seg%thi),Q
       n = ceiling(abs(delta_t/seg%dt))
       delta_t = delta_t/n
       ra = seg%ra(i)
       do i = 1, n
          ra = na_evolve(ra,delta_t)
       end do
       res = one/ra
    end if
    !write(0,'(f15.10,i,f15.12)') t, nseg, res ! HOPPER TESTING
  end function na_Value
  

  !======================================================================
  !! Returns the number of loops in this coupling
  integer function na_NumberOfLoops(nah)
    type(na_handle), intent(in) :: nah
    na_NumberOfLoops = nah%nloop
  end function na_NumberOfLoops
  
  !-------------------------------------------------------------
  !! Do all the necessary cleaning up
  subroutine na_Del(nah)
    type(na_handle), intent(inout) :: nah
    integer :: i, status
    do i = nah%nlo, nah%nhi
       deallocate(nah%seg(i)%ra, stat=status)
    end do
    deallocate(nah%seg, stat=status)
  end subroutine na_Del

  
  !======================================================================
  !! Return the mass of the quark specified by iflv
  real(dp) function na_QuarkMass(nah,iflv) result(mass)
    type(na_handle), intent(in) :: nah
    integer,         intent(in) :: iflv
    if (iflv > 6 .or. iflv < 1) call wae_error('na_QuarkMass', &
         &'illegal value for iflv', intval=iflv)
    mass = nah%quark_masses(iflv)
  end function na_QuarkMass
  

  !-------------------------------------------------------------
  !! Indicate the range of nf values supported by this handle.
  subroutine na_nfRange(nah,nflo,nfhi)
    type(na_handle), intent(in)  :: nah
    integer,         intent(out) :: nflo, nfhi
    if (nah%fixnf == nofixnf) then
       nflo = nah%nlo
       nfhi = nah%nhi
    else
       nflo = nah%fixnf
       nfhi = nah%fixnf
    end if
  end subroutine na_nfRange

  
  !-------------------------------------------------------------
  !! Returns the number of flavours relevant for scale Q
  !! Also returns the range of Qvalues (Qlo, Qhi) for which 
  !! this value of nf remains the same.
  subroutine na_nfAtQ(nah, Q, nfAtQ, Qlo, Qhi, muM_mQ)
    type(na_handle), intent(in)  :: nah
    real(dp),        intent(in)  :: Q
    integer,         intent(out) :: nfAtQ
    real(dp),        intent(out), optional :: Qlo, Qhi
    real(dp),        intent(in),  optional :: muM_mQ
    !-----------------------------
    real(dp) :: deltat_match
    real(dp) :: t
    integer  :: nseg
    real(dp) :: mtlo(nah%nlo:nah%nhi), mthi(nah%nlo:nah%nhi)


    deltat_match = tofQ(default_or_opt(&
         &nah%muMatch_mQuark, muM_mQ)/nah%muMatch_mQuark)

    !-- redefine limits so as to account for requested matching
    !   thresholds muM_mQ. Outer limits should not be modified
    !   since they are rigidly fixed?
    !
    !   What happens if one of these intervals becomes negative?
    mtlo(nah%nlo+1:nah%nhi) = nah%seg(nah%nlo+1:nah%nhi)%tlo + deltat_match
    mthi(nah%nlo:nah%nhi-1) = nah%seg(nah%nlo:nah%nhi-1)%thi + deltat_match
    mtlo(nah%nlo) = nah%seg(nah%nlo)%tlo
    mthi(nah%nhi) = nah%seg(nah%nhi)%thi

    t = tOfQ(Q)
    if (nah%fixnf == nofixnf) then
       do nseg = nah%nlo, nah%nhi
          if (t <= mthi(nseg) .and. t >= mtlo(nseg)) exit
       end do
       if (nseg > nah%nhi) &
            &call wae_Error('na_nfAtQ: Specified Q is not in supported range:'&
            &,dbleval=Q)
       nfAtQ = nseg
       if (present(Qlo) .and. present(Qhi)) then
          Qlo = QOft(mtlo(nseg))
          Qhi = QOft(mthi(nseg))
       end if
    else
       if (t > thi .or. t < tlo) then
          call wae_Error('na_nfAtQ: Specified Q is not in supported range'&
            &,dbleval=Q)
       end if
       nfAtQ = nah%fixnf
       if (present(Qlo) .and. present(Qhi)) then
          Qlo = QOft(tlo)
          Qhi = QOft(thi)
       end if
    end if
  end subroutine na_nfAtQ

  
  !-------------------------------------------------------------
  ! returns the Q range for a given value of nf. If supported.
  subroutine na_QRangeAtNf(nah, nflcl, Qlo, Qhi, muM_mQ)
    type(na_handle), intent(in)  :: nah
    integer,         intent(in)  :: nflcl
    real(dp),        intent(out) :: Qlo, Qhi
    real(dp), optional, intent(in) :: muM_mQ
    !---------------------------------------
    real(dp) :: deltat_match
    character(len=60) :: string

    deltat_match = tofQ(default_or_opt(&
         &nah%muMatch_mQuark, muM_mQ)/nah%muMatch_mQuark)

    if (nah%fixnf == nofixnf) then
       if (nflcl < nah%nlo .or. nflcl > nah%nhi) then
          write(string,'(a,i2,a)') 'nf value ',nflcl,' not supported'
          call wae_Error('QrangeAtNf', trim(string))
       end if
       if (nflcl == nah%nlo) then
          Qlo = QOft(nah%seg(nflcl)%tlo)
       else
          Qlo = QOft(nah%seg(nflcl)%tlo + deltat_match)
       end if
       if (nflcl == nah%nhi) then
          Qhi = QOft(nah%seg(nflcl)%thi)
       else
          Qhi = QOft(nah%seg(nflcl)%thi + deltat_match)
       end if
    else
       if (nflcl /= nah%fixnf) then
          write(string,'(a,i2,a)') 'nf value ',nflcl,' not supported'
          call wae_Error('QrangeAtNf', trim(string))
       end if
       Qlo = QOft(tlo)
       Qhi = QOft(thi)
    end if
  end subroutine na_QRangeAtNf
  
  
  !------------------------------------------------------------
  ! given ra return the value evolve by dt
  function na_evolve(ra,dt) result(res)
    use runge_kutta
    real(dp), intent(in) :: ra, dt
    !----------------------------------
    real(dp) :: res,t
    t = zero
    res = ra
    call rkstp(dt,t,res,na_deriv)
  end function na_evolve
  
  

  !---------------------------------------------------------
  ! derivative of 1/alpha
  subroutine na_deriv(t, ra, dra)
    real(dp), intent(in)  :: t, ra
    real(dp), intent(out) :: dra

    dra = seg%beta0 + seg%beta1/ra + seg%beta2/ra**2
  end subroutine na_deriv

  !-- avoid risk of errors when forgetting factor of two
  function tOfQ(Q)
    real(dp), intent(in) :: Q
    real(dp)             :: tOfQ
    tOfQ = two * log(Q)
  end function tOfQ
  function QOft(t)
    real(dp), intent(in) :: t
    real(dp)             :: QOft
    QOft = exp(half*t)
  end function QOft
  
end module new_as
