!======================================================================
!! Program to run benchmark checks for comparison with Andreas Vogt's
!! numbers
!!
!! Usage: 
!! ./benchmark outname [options]
!!
!! Sends output to files outname.6|7|8|9
!!
!! Options include:
!!
!! Internal choices (accuracy etc.)
!! --------------------------------
!! -dy    dy        spacing in y
!! -dt    dt        4*du
!! -order order     interpolation order
!! -ceps  eps       precision for integrations to get splitting fn.s
!! -lock            whether to lock subgrids together (recommended)
!!
!! -repeat nrep     number of times to repeat calcn (for timing)
!! -precalc         precalculate the evolution
!!
!! -high            write results with high output precision
!! -gnuplot         write output suitable for gnuplot
!!
!! Physics choices
!! ---------------
!! -nloop nloop     # of loops
!! -muR_Q muR_Q     ratio of renorm. to fact. scale
!! -nomasssteps     turn of corrections when going through thresholds
!! -nnlo {exact|param|Nfitav|Nfiterr1|Nfiterr2} 
!!                  which of the NNLO split fns. to use (default param)
!! -nnlo_nfthreshold {exact|param}
!!
!! -varnf           use variable nf
!! -pol             do polarized evln
!!
! $Id: main-tablevogt.f90,v 1.25 2004/09/21 18:50:20 salam Exp $ 
!======================================================================

! the initial conditions
module parton_distributions
  use hoppet_v1
  implicit none
  private
  public :: unpolarized_benchmark_pdf, polarized_benchmark_pdf
contains
  !======================================================================
  !! The dummy PDF suggested by Vogt as the initial condition for the 
  !! unpolazrized evolution
  function unpolarized_benchmark_pdf(xvals) result(pdf)
    real(dp), intent(in) :: xvals(:)
    real(dp)             :: pdf(size(xvals),ncompmin:ncompmax)
    real(dp) :: uv(size(xvals)), dv(size(xvals))
    real(dp) :: ubar(size(xvals)), dbar(size(xvals))
    !---------------------
    real(dp), parameter :: N_g = 1.7_dp, N_ls = 0.387975_dp
    real(dp), parameter :: N_uv=5.107200_dp, N_dv = 3.064320_dp
    real(dp), parameter :: N_db = half*N_ls
  
    pdf = zero
    ! clean method for labelling as PDF as being in the human representation
    call LabelPdfAsHuman(pdf)

    !-- remember that these are all xvals*q(xvals)
    uv = N_uv * xvals**0.8_dp * (1-xvals)**3
    dv = N_dv * xvals**0.8_dp * (1-xvals)**4
    dbar = N_db * xvals**(-0.1_dp) * (1-xvals)**6
    ubar = dbar * (1-xvals)
    pdf(:,iflv_g) = N_g * xvals**(-0.1_dp) * (1-xvals)**5
        
    pdf(:,-iflv_s) = 0.2_dp*(dbar + ubar)
    pdf(:, iflv_s) = pdf(:,-iflv_s)
    pdf(:, iflv_u) = uv + ubar
    pdf(:,-iflv_u) = ubar
    pdf(:, iflv_d) = dv + dbar
    pdf(:,-iflv_d) = dbar
  end function unpolarized_benchmark_pdf


  function polarized_benchmark_pdf(x) result(pdf)
    real(dp), intent(in) :: x(:)
    real(dp)             :: pdf(size(x),-6:7)
    real(dp) :: uv(size(x)), dv(size(x)), ubar(size(x)), dbar(size(x))

    pdf = zero
    ! clean method for labelling as PDF as being in the human representation
    call LabelPdfAsHuman(pdf)
    !-- remember that my definitions that these are all x*q(x)
    uv =  1.3_dp * x**0.7_dp * (1-x)**3 * (1+3*x)
    dv = -0.5_dp * x**0.7_dp * (1-x)**4 * (1+4*x)
    dbar = -0.05_dp * x**0.3_dp * (1-x)**7
    ubar = dbar
    ubar = 0.9_dp * ubar
    dbar = 1.1_dp * dbar
    pdf(:,iflv_g) = 1.5_dp * x**0.5_dp * (1-x)**5

    pdf(:,-iflv_s) = 0.25_dp * (dbar+ubar)
    pdf(:, iflv_s) = pdf(:,-iflv_s)
    pdf(:, iflv_u) = uv + ubar
    pdf(:,-iflv_u) = ubar
    pdf(:, iflv_d) = dv + dbar
    pdf(:,-iflv_d) = dbar

  end function polarized_benchmark_pdf

end module parton_distributions


program main_tablevogt
  ! local things
  use parton_distributions
  use NameSelect; use sub_defs_io
  ! hoppet
  use hoppet_v1
  implicit none
  type(running_coupling)  :: coupling, ash3to4
  integer, parameter :: ngrid = 3
  type(grid_def)     :: grid, gdarray(ngrid)
  real(dp), pointer  :: qdist(:,:)
  type(dglap_holder) :: dh
  real(dp) :: dy, dt, ceps
  integer  :: order, nloop
  logical  :: varnf
  real(dp) :: Qlo
  !------------------------
  !real(dp), parameter :: ymax = 11.6_dp
  real(dp), parameter :: ymax = 18.6_dp
  integer,  parameter :: nms = 20
  real(dp), parameter :: ms(nms) = (/ &
       2.D0,  2.7D0,  3.6D0,   5.D0,   7.D0,   1.D1,&
       1.4D1,  2.D1,   3.D1,   5.D1,   7.D1,   1.D2,&
       2.D2,   5.D2,   1.D3,   3.D3,   1.D4,   4.D4,&
       2.D5,   1.D6 /)
!!$  real(dp), parameter :: ms(20) = (/ &
!!$       2.D0,  2.7D0,  3.6D0,   5.D0,   7.D0,   1.D1,&
!!$       1.4D1,  20.249999999D0,   20.250000001D0,   5.D1,   7.D1,   1.D2,&
!!$       2.D2,   5.D2,   1.D3,   3.D3,   1.D4,   4.D4,&
!!$       2.D5,   1.D6 /)
  real(dp), parameter :: xb(25) = (/ &
       1.D-8,  1.D-7,  1.D-6,&
       1.D-5,  2.D-5,  5.D-5,  1.D-4,  2.D-4,  5.D-4,&
       1.D-3,  2.D-3,  5.D-3,  1.D-2,  2.D-2,  5.D-2,&
       1.D-1, 1.5D-1,  2.D-1,  3.D-1,  4.D-1,  5.0D-1,&
       6.D-1, 7.D-1, 8.D-1, 9.D-1 /)
  type(evln_operator) :: evop(2:nms)

  character(len=*), parameter :: pdfformat_vogt = &
       &'(1X,2(1PE9.1),1X,5(1PE15.7))'
  character(len=*), parameter :: pdfformat_gnu = &
       &'(1X,2(1PE9.1),1X,5(1PE14.6))'
  character(len=*), parameter :: pdfformat_highprec = &
       &'(1X,2(1PE9.1),1X,5(1PE22.14))'
  character(len=len(pdfformat_highprec)) :: pdfformat
  character(len=*), parameter :: asformat = &
       &'(1X,1(1PE9.1),1X,0PF10.7)'
  real(dp) :: qdx(-iflv_max:iflv_max)
  real(dp) :: uv, dv, sv, del, uds, ss, cs, bs, gl, cv
  real(dp) :: x, m2, muR_Q
  integer  :: im, ix, id, nrepeat, irepeat, nfdefault
  integer  :: idev6, idev7, idev8, idev9, idev10
  integer  :: lcl_vogt_imod, factscheme, pdfid
  logical  :: gnuplot, highprec_out, lock, polarized
  logical  :: precalc_evln


  write(0,*) '=============== DISEVEL starting ==================='


  dy    = dble_val_opt('-dy',0.1_dp)
  order = int_val_opt('-order',7)
  nloop = int_val_opt('-nloop',2)
  dt    = dble_val_opt('-dt',0.4_dp)
  ceps  = dble_val_opt('-eps',1e-7_dp)
  muR_Q = dble_val_opt('-muR_Q', one)               !
  muR_Q = sqrt(dble_val_opt('-mu2R_Q', muR_Q**2))   ! allow both specifications
  varnf = log_val_opt('-varnf')
  nfdefault = int_val_opt('-nf',4)
  mass_steps_on = log_val_opt('-masssteps',.true.)
  gnuplot = log_val_opt('-gnuplot')
  highprec_out = log_val_opt('-high')
  lock    = log_val_opt('-lock',.false.)
  nrepeat = int_val_opt('-repeat',1)
  precalc_evln = log_val_opt('-precalc')

  ! take care of consistent use of param and exact variants
  call dglap_Set_nnlo_splitting(&
       &code_val_opt('-nnlo',nnlo_splitting_variant,prefix='nnlo_splitting_'))
  select case(nnlo_splitting_variant)
  case(nnlo_splitting_exact)
     call dglap_Set_nnlo_nfthreshold(nnlo_nfthreshold_exact)
  case default
     call dglap_Set_nnlo_nfthreshold(nnlo_nfthreshold_param)
  end select
  call dglap_Set_nnlo_nfthreshold(&
       &code_val_opt('-nnlo_nfthreshold',nnlo_nfthreshold_variant,&
       &              prefix='nnlo_nfthreshold_'))
  
     

  polarized = log_val_opt('-pol')
  if (polarized) then
     factscheme = factscheme_PolMSbar
  else
     factscheme = factscheme_MSbar
  end if



  call print_comments(prefix='',idev=0)
  idev9 = idev_open_arg(1,'.9')
  !-- security ----------------------
  if (.not. CheckAllArgsUsed(0)) stop
  !----------------------------------
  call print_comments(prefix='# ',idev=idev9)

  if (gnuplot) then
     pdfformat = pdfformat_gnu
  else if (highprec_out) then
     pdfformat = pdfformat_highprec
  else
     pdfformat = pdfformat_vogt
  end if

  call SetDefaultEvolutionDt(dt)
  call SetDefaultConvolutionEps(ceps)

  idev6 = idev_open_arg(1,'.6')
  idev7 = idev_open_arg(1,'.7')
  idev8 = idev_open_arg(1,'.8')
  idev10 = idev_open_arg(1,'.10')


  call InitGridDef(gdarray(3),dy*0.1_dp,0.5_dp, order=order)
  call InitGridDef(gdarray(2),dy/3.0_dp,2.0_dp, order=order)
  call InitGridDef(gdarray(1),dy,       ymax  ,order=order)
  call InitGridDef(grid,gdarray(1:3),locked=lock)
  
  call qcd_SetNf(nfdefault) ! <-- CHECK THIS
  if (varnf) then
     ! given that we do not actually need to evolve in nf=3, but just
     ! go past the threshold, we could not bother with this nf=3,
     ! but it makes things a bit simpler.
     call InitDglapHolder(grid, dh, factscheme, nloop=nloop,&
          & nflo=3, nfhi=6)
     if (muR_Q < one) then
        !-- this is a hack to make sure that alfas=0.35 corresponds
        !   to the 3-flavour alphas.
        !   initialise with 3 flavours (default for this scale)
        call InitRunningCoupling(ash3to4, alfas=0.35_dp, Q=sqrt(ms(1)), nloop=nloop,&
             & muMatch_mQuark=one, use_nah=.true.)
        !   reinitialise with new alpha_s with muMatch_mQuark=muR_Q
        !   and 4 flavours
        !call InitRunningCoupling(coupling, alfas=Value(ash3to4,Q=sqrt(ms(1)),fixnf=4),&
        !     & Q=sqrt(ms(1)), nloop=nloop,&
        !     & muMatch_mQuark=muR_Q, use_nah=.true.)
        !------------------
        ! "second scenario:" evolve down to muR_Q*mc and set the
        ! new coupling there.
        call InitRunningCoupling(coupling, alfas=Value(ash3to4,Q=sqrt(ms(1))*muR_Q),&
             & Q=sqrt(ms(1))*muR_Q, nloop=nloop,&
             & muMatch_mQuark=muR_Q, use_nah=.true.)
     else
        call InitRunningCoupling(coupling, alfas=0.35_dp, Q=sqrt(ms(1)), nloop=nloop,&
             & muMatch_mQuark=muR_Q, use_nah=.true.)
     end if
     
!!$     call InitRunningCoupling(coupling, alfas=0.115161, Q=100.0_dp, nloop=nloop,&
!!$          & use_nah=.true.)
  else
     call InitDglapHolder(grid, dh, factscheme, nloop=nloop)
     call InitRunningCoupling(coupling, alfas=0.35_dp, Q=sqrt(ms(1)), nloop=nloop,&
          & fixnf=nf_int, use_nah=.true.)
  end if
  

  write(0,*) 'Init done'

  Qlo = sqrt(ms(1))

  !--- allow testing of precalculated and stored evolution
  if (precalc_evln) then
     do im = 2, size(ms)
        call EvolveGeneric(dh, nloop=nloop, coupling=coupling, &
             &Q_init=sqrt(ms(im-1)), Q_end=sqrt(ms(im)), &
             &evop = evop(im), muR_Q = muR_Q)
     end do
  end if

  do irepeat = 1, nrepeat
     
  call allocPDF(grid, qdist)
  if (polarized) then
     qdist = polarized_benchmark_pdf(xValues(grid))
  else 
     qdist = unpolarized_benchmark_pdf(xValues(grid))
  end if

  do im = 1, size(ms)
     m2 = ms(im)
     if (im > 1) then
        if (precalc_evln) then
           qdist = evop(im) .conv. qdist
        else
           call EvolvePDF(dh, qdist, nloop=nloop, coupling=coupling, &
                &Q_init=sqrt(ms(im-1)), Q_end=sqrt(ms(im)), muR_Q = muR_Q)
        end if
     end if
     
     !if (irepeat /= 1) cycle
     if (irepeat == 1) write(idev6,asformat) m2, Value(coupling, sqrt(m2))
     if (irepeat > 1) cycle ! get JUST the time for evolution...
                            ! later on worry about optimizing the rest...
     do ix = 1, size(xb)
        x = xb(ix)
        !-- decant...
        qdx(iflv_min:iflv_max) = &
             &EvalGridQuant(grid,qdist(:,iflv_min:iflv_max),-log(x))
        uv  = qdx(iflv_u) - qdx(iflv_ubar)
        dv  = qdx(iflv_d) - qdx(iflv_dbar)
        sv  = qdx(iflv_s) - qdx(iflv_sbar); 
        del = qdx(iflv_dbar) - qdx(iflv_ubar)
        uds = two*(qdx(iflv_dbar) + qdx(iflv_ubar)) ! presume this is meant
        ss  = qdx(iflv_s) + qdx(iflv_sbar)
        cs  = qdx(iflv_c) + qdx(iflv_cbar)
        bs  = qdx(iflv_b) + qdx(iflv_bbar)
        gl  = qdx(iflv_g)
        cv  = qdx(iflv_c) - qdx(iflv_cbar); !TEMPORARY FOR DEBUGGING
        !-- fix up rounding errors to avoid embarassment!
        if (nloop <= 2) sv = zero
        if (irepeat == 1) then
           write(idev7,pdfformat) M2, X, UV, DV, DEL, UDS, GL
           write(idev8,pdfformat) M2, X, SV, SS, CS, BS
           if (abs(m2-1e4_dp) < 1e-3_dp) call vogt_latex
        end if
     end do
     if (gnuplot .and. irepeat == 1) then
        write(idev7,*) 
        write(idev8,*) 
        write(idev7,*) 
        write(idev8,*) 
     end if
  end do
  deallocate(qdist)
  end do


contains

  subroutine print_comments(prefix,idev)
    character(len=*), intent(in) :: prefix
    integer,          intent(in) :: idev

    character(len=*), parameter :: fa='(a,a)'
    character(len=*), parameter :: fd='(a,f12.8)'
    character(len=*), parameter :: fe='(a,es15.5)'
    character(len=*), parameter :: fi='(a,i6)'
    character(len=*), parameter :: fl='(a,l2)'

    write(idev,fa) prefix//'Output will be sent to '//trim(string_val_arg(1))//'.[678]'
    write(idev,fd) prefix//'dy = ', dy
    write(idev,fd) prefix//'dt = ', dt
    write(idev,fi) prefix//'num order = ', order
    write(idev,fl) prefix//'locking is ', lock
    write(idev,fi) prefix//'nloop = ', nloop
    write(idev,fe) prefix//'conv_eps = ', ceps
    write(idev,fd) prefix//'muR_Q = ', muR_Q
    write(idev,fd) prefix//'muR_Q^2 = ', muR_Q**2
    write(idev,fl) prefix//'varnf is ',varnf
    if (.not. varnf) &
         & write(idev,fi) prefix//'nf = ',nfdefault
    write(idev,fl) prefix//'polarization = ', polarized
    write(idev,fl) prefix//'mass_steps_on is',mass_steps_on
    write(idev,fa) prefix//'factscheme = ', &
         &trim(NameOfCode(factscheme,prefix='factscheme'))
    if (nloop > 2) then
       write(idev,fa) prefix//'nnlo_splitting_variant = ', &
         &trim(NameOfCode(nnlo_splitting_variant,prefix='nnlo_splitting'))
       write(idev,fa) prefix//'nnlo_nfthreshold_variant = ', &
         &trim(NameOfCode(nnlo_nfthreshold_variant,prefix='nnlo_nfthreshold'))
    end if
    write(idev,fl) prefix//'Precalc evolution   ', precalc_evln
    write(idev,fi) prefix//'Number of repeats ', nrepeat
    write(idev,fa,advance="no") prefix
    call time_stamp(idev)
    write(idev,fa) '# '//trim(command_line())
  end subroutine print_comments
  

  subroutine vogt_latex
    !DATA XB / 1.D-8,  1.D-7,  1.D-6,  1.D-5,  1.D-4,  1.D-3,         &
    !     &           1.D-2,  1.D-1,  3.D-1,  5.D-1,  7.D-1,  9.D-1 /        
    !                                                                       
    ! ..Output                                                              
    !                                                                       
    !        WRITE (7,18) M2, X, UV, DV, DEL, UDS, GL                       
    !        WRITE (8,18) M2, X, SV, SS, CS, BS                             
    !WRITE (7,18) X, UV, DV, DEL, UDS, GL 
    !WRITE (8,18) X, SV, SS, CS, BS 
    !                                                                       
    logical :: writeit

    writeit = .false.
    IF ( (X .GT. 0.9E-7) .and. (X .LT. 1.1E-7) ) THEN
       writeit = .true. ; WRITE (idev10,20)
    ELSE IF ( (X .GT. 0.9E-6) .and. (X .LT. 1.1E-6) ) THEN
       writeit = .true. ; WRITE (idev10,21)
    ELSE IF ( (X .GT. 0.9E-5) .and. (X .LT. 1.1E-5) ) THEN
       writeit = .true. ; WRITE (idev10,22)
    ELSE IF ( (X .GT. 0.9E-4) .and. (X .LT. 1.1E-4) ) THEN
       writeit = .true. ; WRITE (idev10,23)
    ELSE IF ( (X .GT. 0.9E-3) .and. (X .LT. 1.1E-3) ) THEN
       writeit = .true. ; WRITE (idev10,24)
    ELSE IF ( (X .GT. 0.9E-2) .and. (X .LT. 1.1E-2) ) THEN
       writeit = .true. ; WRITE (idev10,25)
    ELSE IF ( (X .GT. 0.9E-1) .and. (X .LT. 1.1E-1) ) THEN
       writeit = .true. ; WRITE (idev10,26)
    ELSE IF ( (X .GT. 2.9E-1) .and. (X .LT. 3.1E-1) ) THEN
       writeit = .true. ; WRITE (idev10,27)
    ELSE IF ( (X .GT. 4.9E-1) .and. (X .LT. 5.1E-1) ) THEN
       writeit = .true. ; WRITE (idev10,28)
    ELSE IF ( (X .GT. 6.9E-1) .and. (X .LT. 7.1E-1) ) THEN
       writeit = .true. ; WRITE (idev10,29)
    ELSE IF ( (X .GT. 8.9E-1) .and. (X .LT. 9.1E-1) ) THEN
       writeit = .true. ; WRITE (idev10,30)
    END IF
    if (writeit) then
       WRITE (idev10,31) UV, DV, DEL, UDS
       ! WRITE (idev10,32) SV, SS, CS, GL
       WRITE (idev10,32) SS, CS, BS, GL
    end if
    
20  FORMAT ('$q! 10^{-7}q!$ &')
21  FORMAT ('$q! 10^{-6}q!$ &')
22  FORMAT ('$q! 10^{-5}q!$ &')
23  FORMAT ('$q! 10^{-4}q!$ &')
24  FORMAT ('$q! 10^{-3}q!$ &')
25  FORMAT ('$q! 10^{-2}q!$ &')
26  FORMAT ('$q! 0.1    q!$ &')
27  FORMAT ('$q! 0.3    q!$ &')
28  FORMAT ('$q! 0.5    q!$ &')
29  FORMAT ('$q! 0.7    q!$ &')
30  FORMAT ('$q! 0.9    q!$ &')
31  FORMAT (3('$q!',1PE11.4,'}q!$ &'),'$q!',1PE11.4,'}q!$ &')
32  FORMAT (3('$q!',1PE11.4,'}q!$ &'),'$q!',1PE11.4,'}q!$qq')

!!$    IF ( (X .GT. 0.9E-7) .and. (X .LT. 1.1E-7) ) THEN 
!!$       WRITE (idev10,20) UV, DV, DEL, UDS 
!!$       WRITE (idev10,40) SS, CS, GL 
!!$    ELSE IF ( (X .GT. 0.9E-6) .and. (X .LT. 1.1E-6) ) THEN 
!!$       WRITE (idev10,21) UV, DV, DEL, UDS 
!!$       WRITE (idev10,41) SS, CS, GL 
!!$    ELSE IF ( (X .GT. 0.9E-5) .and. (X .LT. 1.1E-5) ) THEN 
!!$       WRITE (idev10,22) UV, DV, DEL, UDS 
!!$       WRITE (idev10,42) SS, CS, GL 
!!$    ELSE IF ( (X .GT. 0.9E-4) .and. (X .LT. 1.1E-4) ) THEN 
!!$       WRITE (idev10,23) UV, DV, DEL, UDS 
!!$       WRITE (idev10,43) SS, CS, GL 
!!$    ELSE IF ( (X .GT. 0.9E-3) .and. (X .LT. 1.1E-3) ) THEN 
!!$       WRITE (idev10,24) UV, DV, DEL, UDS 
!!$       WRITE (idev10,44) SS, CS, GL 
!!$    ELSE IF ( (X .GT. 0.9E-2) .and. (X .LT. 1.1E-2) ) THEN 
!!$       WRITE (idev10,25) UV, DV, DEL, UDS 
!!$       WRITE (idev10,45) SS, CS, GL 
!!$    ELSE IF ( (X .GT. 0.9E-1) .and. (X .LT. 1.1E-1) ) THEN 
!!$       WRITE (idev10,26) UV, DV, DEL, UDS 
!!$       WRITE (idev10,46) SS, CS, GL 
!!$    ELSE IF ( (X .GT. 2.9E-1) .and. (X .LT. 3.1E-1) ) THEN 
!!$       WRITE (idev10,27) UV, DV, DEL, UDS 
!!$       WRITE (idev10,47) SS, CS, GL 
!!$    ELSE IF ( (X .GT. 4.9E-1) .and. (X .LT. 5.1E-1) ) THEN 
!!$       WRITE (idev10,28) UV, DV, DEL, UDS 
!!$       WRITE (idev10,48) SS, CS, GL 
!!$    ELSE IF ( (X .GT. 6.9E-1) .and. (X .LT. 7.1E-1) ) THEN 
!!$       WRITE (idev10,29) UV, DV, DEL, UDS 
!!$       WRITE (idev10,49) SS, CS, GL 
!!$    ELSE IF ( (X .GT. 8.9E-1) .and. (X .LT. 9.1E-1) ) THEN 
!!$       WRITE (idev10,30) UV, DV, DEL, UDS 
!!$       WRITE (idev10,50) SS, CS, GL 
!!$    END IF
!!$    !                                                                       
!!$10  FORMAT (1PE8.1,1X,7(1PE10.3)) 
!!$12  FORMAT (1X,1PE9.1,1X,6(1PE11.4)) 
!!$18  FORMAT (1X,1(1PE9.1),1X,5(1PE13.6)) 
!!$    !                                                                       
!!$20  FORMAT (1X,'$10^{-7}$ &',3(1PE11.4,'}$ &'),1PE11.4,'}$ ') 
!!$40  FORMAT (1X,'          &',2(1PE11.4,'}$ &'),1PE11.4,'}$ qq') 
!!$21  FORMAT (1X,'$10^{-6}$ &',3(1PE11.4,'}$ &'),1PE11.4,'}$ ') 
!!$41  FORMAT (1X,'          &',2(1PE11.4,'}$ &'),1PE11.4,'}$ qq') 
!!$22  FORMAT (1X,'$10^{-5}$ &',3(1PE11.4,'}$ &'),1PE11.4,'}$ ') 
!!$42  FORMAT (1X,'          &',2(1PE11.4,'}$ &'),1PE11.4,'}$ qq') 
!!$23  FORMAT (1X,'$10^{-4}$ &',3(1PE11.4,'}$ &'),1PE11.4,'}$ ') 
!!$43  FORMAT (1X,'          &',2(1PE11.4,'}$ &'),1PE11.4,'}$ qq') 
!!$24  FORMAT (1X,'$10^{-3}$ &',3(1PE11.4,'}$ &'),1PE11.4,'}$ ') 
!!$44  FORMAT (1X,'          &',2(1PE11.4,'}$ &'),1PE11.4,'}$ qq') 
!!$25  FORMAT (1X,'$10^{-2}$ &',3(1PE11.4,'}$ &'),1PE11.4,'}$ ') 
!!$45  FORMAT (1X,'          &',2(1PE11.4,'}$ &'),1PE11.4,'}$ qq') 
!!$26  FORMAT (1X,'0.1       &',3(1PE11.4,'}$ &'),1PE11.4,'}$ ') 
!!$46  FORMAT (1X,'          &',2(1PE11.4,'}$ &'),1PE11.4,'}$ qq') 
!!$27  FORMAT (1X,'0.3       &',3(1PE11.4,'}$ &'),1PE11.4,'}$ ') 
!!$47  FORMAT (1X,'          &',2(1PE11.4,'}$ &'),1PE11.4,'}$ qq') 
!!$28  FORMAT (1X,'0.5       &',3(1PE11.4,'}$ &'),1PE11.4,'}$ ') 
!!$48  FORMAT (1X,'          &',2(1PE11.4,'}$ &'),1PE11.4,'}$ qq') 
!!$29  FORMAT (1X,'0.7       &',3(1PE11.4,'}$ &'),1PE11.4,'}$ ') 
!!$49  FORMAT (1X,'          &',2(1PE11.4,'}$ &'),1PE11.4,'}$ qq') 
!!$30  FORMAT (1X,'0.9       &',3(1PE11.4,'}$ &'),1PE11.4,'}$ ') 
!!$50  FORMAT (1X,'          &',2(1PE11.4,'}$ &'),1PE11.4,'}$ qq')
  end subroutine vogt_latex
  
end program main_tablevogt
