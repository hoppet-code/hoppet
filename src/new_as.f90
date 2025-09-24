module new_as
  use types; use consts_dp
  use assertions; use warnings_and_errors
  implicit none
  private

  type na_segment
     real(dp) :: tlo, thi, dt, invdt
     real(dp) :: beta0, beta1, beta2, beta3
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
     logical                   :: masses_are_MSbar
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
  public :: na_Init, na_Value, na_Value_full, na_Del, na_NumberOfLoops
  public :: na_nfRange, na_nfAtQ, na_QrangeAtNf, na_QuarkMass
  public :: na_QuarkMassesAreMSbar
  public :: na_Set_dt_base

  !type(na_segment), pointer :: seg

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
  
  subroutine na_Init(nah, alfas, Q, nloop, fixnf, quark_masses, masses_are_MSbar, muMatch_mQuark)
    use qcd; 
    type(na_handle),  intent(out), target :: nah
    real(dp), intent(in), optional :: alfas, Q
    integer,  intent(in), optional :: nloop, fixnf
    real(dp), intent(in), optional :: quark_masses(4:6)
    real(dp), intent(in), optional :: muMatch_mQuark
    logical,  intent(in), optional :: masses_are_MSbar
    !---------------------------------------------
    !real(dp) :: alfas_lcl, Q_lcl
    !integer  :: nloop_lcl, fixnf_lcl
    !-------------------------------
    real(dp) :: alfas_here, lnmatch, alphastep20_lclscheme
    real(dp) :: alphastep30_lclscheme, alphastep31_lclscheme 
    real(dp) :: alphastep31_inv_lclscheme 
    real(dp) :: alphastep30_nl_lclscheme, alphastep31_nl_lclscheme 
    real(dp) :: tstart, t, ra
    integer  :: nf_store
    integer  :: nbin, i, j, nseg, istart
    integer  :: nbin_total
    type(na_segment), pointer :: seg

    !-- we may well play with nf, so need to be able to reset it
    nf_store = nf_int

    nah%alfas = default_or_opt(0.118_dp, alfas)
    nah%Q     = default_or_opt(91.2_dp, Q)
    nah%nloop = default_or_opt(2, nloop)

    nah%fixnf = default_or_opt(nofixnf, fixnf)

    nah%muMatch_mQuark = default_or_opt(one, muMatch_mQuark)
    nah%masses_are_MSbar = default_or_opt(.false., masses_are_MSbar)

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

    ! check we have CA=3 and CF=4/3 for the nloop >=4 case
    if (nah%nloop >= 4) then
       if (ca /= ca_def) then
          call wae_error('na_Init','CA/=CA_def not supported for nloop >=4. CA is',dbleval=ca)
       else if (abs(cf - cf_def) > 1e-10_dp) then
          call wae_error('na_Init','CF/=CF_def not supported for nloop >=4. CF is',dbleval=cf)
       else if (abs(tr - tr_def) > 1e-10_dp) then
          call wae_error('na_Init','TR/=TR_def not supported for nloop >=4. TR is',dbleval=tr)
       end if
    end if
    
    
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
       case(0)
         seg%beta0 = zero
         seg%beta1 = zero
         seg%beta2 = zero
         seg%beta3 = zero
       case(1)
          seg%beta1 = zero
          seg%beta2 = zero
          seg%beta3 = zero
       case(2)
          seg%beta1 = beta1
          seg%beta2 = zero
          seg%beta3 = zero
       case(3)
          seg%beta1 = beta1
          seg%beta2 = beta2
          seg%beta3 = zero
       case(4)
          seg%beta1 = beta1
          seg%beta2 = beta2
          seg%beta3 = beta3
       case default
          call wae_Error('na_init: unsupported number of loops requested')
       end select

       !-- first guess, which then needs to be "adapted" to the range
       seg%dt    = dt_base * half**(6-i)
       nbin = ceiling((seg%thi - seg%tlo)/seg%dt)
       nbin_total = nbin_total + nbin
       seg%dt    = (seg%thi - seg%tlo) / nbin
       seg%invdt = one / seg%dt
       allocate(seg%ra(0:nbin))
    end do
    !write(0,*) 'new_as: used total bins ', nbin_total


    ! nloop = 0 will be fixed coupling, in which case there is nothing
    ! else to do (we put this after the last of the "allocate"
    ! statements, to help ensure that deallocation works smoothly)
    if (nah%nloop == 0) return

    !-- find out in which segment we start
    tstart = tOfQ(nah%Q)
    do i = nah%nlo, nah%nhi
       if (tstart <= nah%seg(i)%thi .and. tstart >= nah%seg(i)%tlo) exit
    end do
    if (i > nah%nhi) &
         &call wae_Error('na_init: Specified Q not in supported range')
    nseg = i

    !-- fill up the starting segment
    ra = one/nah%alfas
    seg => nah%seg(nseg)
    istart = nint((tstart-seg%tlo)/seg%dt)
    t = tstart
    seg%ra(istart) = na_evolve(ra, seg%tlo+istart*seg%dt - tstart,seg)
    do i = istart+1, ubound(seg%ra,dim=1)
       seg%ra(i) = na_evolve(seg%ra(i-1), seg%dt, seg)
    end do
    do i = istart-1, lbound(seg%ra,dim=1), -1
       seg%ra(i) = na_evolve(seg%ra(i+1), -seg%dt, seg)
    end do
    seg%iflip = istart

    ! handle pole v. MSbar mass treatment
    if (nah%masses_are_MSbar) then
      alphastep20_lclscheme = alphastep20_msbar
      alphastep30_lclscheme = alphastep30_msbar
      alphastep31_lclscheme = alphastep31_msbar
      alphastep31_inv_lclscheme = alphastep31_msbar_inv
      alphastep30_nl_lclscheme = alphastep30_msbar_nl
      alphastep31_nl_lclscheme = alphastep31_msbar_nl
    else
      alphastep20_lclscheme = alphastep20_pole
      alphastep30_lclscheme = alphastep30_pole
      alphastep31_lclscheme = alphastep31_pole
      alphastep31_inv_lclscheme = alphastep31_pole_inv
      alphastep30_nl_lclscheme = alphastep30_pole_nl
      alphastep31_nl_lclscheme = alphastep31_pole_nl
    end if

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
                  &    + alphastep20_lclscheme)*alfas_here**2)
          case(4)
             alfas_here = alfas_here*(1 + alphastep11 * lnmatch * alfas_here&
                  & + (alphastep22 * lnmatch**2 + alphastep21 * lnmatch&
                  & + alphastep20_lclscheme)*alfas_here**2&
                  & + (alphastep33 * lnmatch**3 + alphastep32_inv * lnmatch**2&
                  & + (alphastep31_inv_lclscheme - (j-1) * alphastep31_nl_lclscheme) * lnmatch & 
                  & + (-alphastep30_lclscheme - (j-1) * alphastep30_nl_lclscheme))*alfas_here**3)
          end select
       end if
       seg%ra(0) = one/alfas_here
       ! recall that this is the reciprocal of alpha!
       do i = 1, ubound(seg%ra,dim=1)
          seg%ra(i) = na_evolve(seg%ra(i-1), seg%dt, seg)
       end do
       seg%iflip = -1
    end do
    
    !-- fill up the segments below
    do j = nseg-1, nah%nlo, -1
       seg => nah%seg(j)
       alfas_here = one/nah%seg(j+1)%ra(0)
       if (mass_steps_on .and. nah%fixnf == nofixnf) then
!          write(0,*) '------ changing as from ',j+1,' to ',j
!          write(0,*) j+1, alfas_here
          lnmatch = two*log(nah%muMatch_mQuark)
          select case (nah%nloop)
          case(2)
             alfas_here = alfas_here*(1 - alphastep11 * lnmatch * alfas_here)
          case(3)
             alfas_here = alfas_here*(1 - alphastep11 * lnmatch * alfas_here&
                  & + (alphastep22 * lnmatch**2 - alphastep21 * lnmatch&
                  &    - alphastep20_lclscheme)*alfas_here**2)
          case(4)
             alfas_here = alfas_here*(1 - alphastep11 * lnmatch * alfas_here&
                  & + (alphastep22 * lnmatch**2 - alphastep21 * lnmatch&
                  & - alphastep20_lclscheme)*alfas_here**2&
                  & + (-alphastep33 * lnmatch**3 - alphastep32 * lnmatch**2&
                  & + (-alphastep31_lclscheme + j * alphastep31_nl_lclscheme) * lnmatch & 
                  & + (alphastep30_lclscheme + j * alphastep30_nl_lclscheme))*alfas_here**3)
          end select
!          write(0,*) j, alfas_here
!          write(0,*) 'Coefficients', alphastep20_lclscheme * pisq
       end if
       seg%ra(ubound(seg%ra,dim=1)) = one/alfas_here
       do i = ubound(seg%ra,dim=1)-1, 0, -1
          seg%ra(i) = na_evolve(seg%ra(i+1), -seg%dt, seg)
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
  function na_Value_full(nah, Q, fixnf) result(res)
    type(na_handle), intent(in), target :: nah
    real(dp),        intent(in)         :: Q
    integer,         intent(in), optional :: fixnf
    real(dp)                            :: res
    !-------------------------------------
    real(dp) :: t, trescaled, ra, delta_t
    integer  :: nseg, i, n
    integer, save :: warn_id = warn_id_INIT
    integer, parameter :: max_warn = 1
    type(na_segment), pointer :: seg

    if (nah%nloop == 0) then
      res = nah%alfas
      return
    end if

    t = tOfQ(Q)

    if (present(fixnf) .and. nah%fixnf == nofixnf) then
       if (fixnf < nah%nlo .or. fixnf > nah%nhi) then
          call wae_error('na_Value_full:', 'the fixnf requested is&
               & outside the range supported this na_handle')
       end if
       if (t < tlo .or. t > thi) then
          call wae_error('na_Value_full:', 'the Q value is&
               & outside the range supported this na_handle')
       end if
       nseg = fixnf
       !if (t > nah%seg(nseg)%thi+ nah%seg(nseg)%dt .or.&
       !     & t < nah%seg(nseg)%tlo-nah%seg(nseg)%dt) then
       !   call wae_error('na_Value_full:', &
       !        &  'With fixnf, Q is too far outside supported range.')
       !end if
    else
       if (present(fixnf) .and. nah%fixnf /= nofixnf) then
          if (fixnf /= nah%fixnf) then
             call wae_error('na_Value_full:', 'the fixnf requested is &
                  &different from that supported by na_handle')
          end if
       end if
       do nseg = nah%nlo, nah%nhi
          if (t <= nah%seg(nseg)%thi .and. t >= nah%seg(nseg)%tlo) exit
       end do
       if (nseg > nah%nhi) &
            &call wae_Error('na_Value_full: Specified Q is not in supported range'&
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
       res = one/na_evolve(seg%ra(i), delta_t, seg)
    else
       call wae_warn(max_warn,warn_id,'na_Value_full: will evolve &
            &fixed-nf alpha_s beyond precalculated range.',&
            &'This procedure may be very slow, Q=', dbleval=Q)
       !write(0,*) Qoft(seg%tlo),Qoft(seg%thi),Q
       n = ceiling(abs(delta_t/seg%dt))
       delta_t = delta_t/n
       ra = seg%ra(i)
       do i = 1, n
          ra = na_evolve(ra,delta_t, seg)
       end do
       res = one/ra
    end if
    !write(0,'(f15.10,i,f15.12)') t, nseg, res ! HOPPER TESTING
  end function na_Value_full

  !! Faster evaluation of alpha_s, with fewer options; reverts
  !! to _full for some edge cases and also doesn't have the 
  !! fixnf option
  function na_Value(nah, Q, fixnf) result(res)
    use interpolation_coeffs; use interpolation
    type(na_handle), intent(in), target :: nah
    real(dp),        intent(in)         :: Q
    integer,         intent(in), optional :: fixnf
    real(dp) :: res
    !---------------
    real(dp) :: t, tnorm, tdarr(0:4), itd, prod
    ! on M2Pro-gfortran15-O3 there's a small (~2%, 0.2ns out of 10ns) speed advantage
    ! in having the interpolation coefficients hard-coded here
    real(dp), parameter :: coeffs(0:4) = [one/24.0_dp, -one/6.0_dp, one/4.0_dp, -one/6.0_dp, one/24.0_dp]
    integer  :: iseg, it, i
    type(na_segment), pointer :: this_seg
    
    if (present(fixnf)) then
      res = na_Value_full(nah, Q, fixnf)
      return
    end if

    t = tofQ(Q)
    if (t < nah%seg(nah%nlo)%tlo .or. t > nah%seg(nah%nhi)%thi) then
      res = na_Value_full(nah,Q)
      return
    end if

    ! choose the segment
    do iseg = nah%nlo, nah%nhi
      if (t <= nah%seg(iseg)%thi) exit
    end do

    this_seg => nah%seg(iseg)
    if (ubound(this_seg%ra,1) < 4) then
      ! not enough points to do the interpolation, so just use na_value_full
      res = na_Value_full(nah,Q)
      return
    end if

    ! 2025-09-06: on M2Pro -O3 gfortran, save about 1ns by using multiplication rather than division
    tnorm = (t - this_seg%tlo) * this_seg%invdt
    it = int(tnorm)
    if (it == tnorm) then
      res = one/this_seg%ra(it)
      return
    end if

    ! M2Pro gfortran 15 -O3, save about 0.5 ns with explicit if rather than min,max
    ! So do not use `it = min(max(it - 2,0),ubound(this_seg%ra,1)-4)`
    if (it < 0) then
      it = it + 2
    else if (it + 4 > ubound(this_seg%ra,1)) then
      it = ubound(this_seg%ra,1) - 4
    end if

    itd = it
    ! This route gives a total timing of 10.0-10.2 ns on M2Pro-gfortran15-O3
    tdarr = (tnorm - itd) - [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    !!! this route involves a division, but on M2Pro-gfortran15-O3 save 0.2-0.3ns
    prod = product(tdarr)
    res = prod * sum((coeffs * this_seg%ra(it:it+4)) / tdarr)

    ! this route skips divisions but on M2Pro-gfortran15-O3 costs an extra 0.2-0.3ns
    !block
    !  real(dp) :: prodsl(0:4), prodsr(0:4)
    !  prodsl(0) = one
    !  prodsr(4) = one
    !  do i = 1,4
    !    prodsl(i) = prodsl(i-1) * tdarr(i-1)
    !    prodsr(4-i) = prodsr(5-i) * tdarr(5-i)
    !  end do
    !  res = sum(prodsl * prodsr * this_seg%ra(it:it+4) * coeffs)
    !end block 

    !! alt1: 12.2-12.4 ns
    !call fill_interp_weights4(tnorm-itd, tdarr)
    !! alt2: 14.7-15.0 ns
    !call uniform_interpolation_weights(tnorm-itd, tdarr)
    !res = sum(tdarr * this_seg%ra(it:it+4))

    res = one / res
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
  ! return true if the quark masses are in the MSbar scheme
  logical function na_QuarkMassesAreMSbar(nah) result(res)
    type(na_handle), intent(in) :: nah
    res = nah%masses_are_MSbar
  end function na_QuarkMassesAreMSbar
  

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
  function na_evolve(ra,dt, seg) result(res)
    use runge_kutta
    real(dp), intent(in) :: ra, dt
    class(*), intent(in) :: seg
    !----------------------------------
    real(dp) :: res,t
    t = zero
    res = ra
    !call rkstp_arg(dt,t,res,na_deriv,seg)
    call rkstp_arg(dt,t,res,na_deriv,seg)
  end function na_evolve
  
  !---------------------------------------------------------
  !! derivative of 1/alpha
  subroutine na_deriv(t, ra, dra, seg)
    real(dp), intent(in)  :: t, ra
    real(dp), intent(out) :: dra
    class(*), intent(in)  :: seg

    select type(seg)
    type is (na_segment)
      dra = seg%beta0 + seg%beta1/ra + seg%beta2/ra**2 + seg%beta3/ra**3
    class default
      call wae_error('na_deriv: unknown polymorphic type for seg')
    end select
  end subroutine na_deriv

  !! Function that returns t = 2*log(Q).
  !! It is coded as a separate function to avoid the risk of forgetting
  !! the factor of 2 in the log and to allow for checks that Q is not zero.
  function tOfQ(Q)
    real(dp), intent(in) :: Q
    real(dp)             :: tOfQ
    if (Q <= zero) then
      ! anything that is too small is set to a very negative value for tofQ;
      ! Recall: huge(Q) is the largest representable number of the same type as Q.
      tOfQ = -huge(Q)
    else
      tOfQ = two * log(Q)
    end if
  end function tOfQ

  !! Function that returns Q = exp(half*t).
  function QOft(t)
    real(dp), intent(in) :: t
    real(dp)             :: QOft
    QOft = exp(half*t)
  end function QOft
  
end module new_as
