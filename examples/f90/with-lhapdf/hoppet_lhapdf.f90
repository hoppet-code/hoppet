
!! Module to facilitate the transfer of LHAPDF PDF sets into HOPPET
!!
!! The main routines are:
!!
!! - LoadLHAPDF(name[,imem=0]): 
!!   Load a single LHAPDF member into streamlined interface (tables(0))
!!
!! - LoadLHAPDFSet(name):
!!   Load all members of an LHAPDF set, with member 0 going into tables(0)
!!   of the streamlined interface, and all members going into
!!   lhapdf_tables(0:max_mem) (where max_mem is the maximum member number
!!   of the set).  
!!
!!
module hoppet_lhapdf
  use streamlined_interface, InitPDF_hoppet => InitPDF
  use hoppet_term
  use hoppet_lhapdf_interfaces
  use types
  use pdf_tabulate
  use qcd_coupling
  use hoppet_to_string
  implicit none

  private

  !! a type to hold information about an LHAPDF set loaded into HOPPET
  type hoppet_lhapdf_set
    !! range of the original LHAPDF set
    real(dp) :: Qmin, Qmax, xmin, xmax

    !! quark masses for VFN scheme
    real(dp) :: mc, mb, mt

    !! index of the highest available member in the set
    integer  :: max_mem = -1

    !! name of the LHAPDF set (potentially truncated)
    character(len=300) :: name = ""

    !! PDF set metadata
    integer  :: orderPDF  !< orderPDF=1 means LO evolution
    integer  :: nloop     !< number of loops from a hoppet point of view
    integer  :: nfmax     !< maximum number of flavours in the set

    type(pdf_table),        pointer :: tables(:)    => null() ! array of tables
    type(running_coupling), pointer :: couplings(:) => null() ! array of couplings
    logical :: owns_data = .false.
  end type hoppet_lhapdf_set

  type(hoppet_lhapdf_set), target, save :: lhapdf_set

  character(len=*), parameter :: tfmt = '("'//blue//'",a,f8.2,a,"'//reset//'")'


  public :: LoadLHAPDF, LoadLHAPDFSet
  public :: hoppet_lhapdf_set
  public :: lhapdf_set

contains

  !! Load a single LHAPDF member into streamlined interface (tables(0))
  !! If imem is not given, member 0 is loaded
  !!
  !! The routine also starts up hoppet with the right range and sensible
  !! default parameters
  !! 
  subroutine LoadLHAPDF(name, imem)
    character(len=*),  intent(in) :: name
    integer, optional, intent(in) :: imem
    !---------------------------------------
    real(dp) :: Q2minPDF, Q2maxPDF
    type(hoppet_lhapdf_set), pointer :: set
    real(dp) :: ta1, ta2
    real(dp) :: ymax, dy, dlnlnQ!, Q0
    integer  :: order, yorder, lnlnQorder
    real(dp), parameter :: mz = 91.1880_dp

    set => lhapdf_set ! shorthand

    ! Load LHAPDF set
    call cpu_time(ta1)
    call initPDFSetByName(trim(name))
    if (present(imem)) then
      call initPDF(imem)
    else
      call initPDF(0)
    end if

    call numberpdf(set%max_mem)
    !call cpu_time(ta2)
    !write(*,tfmt) "Time to load LHAPDF sets: ", (ta2-ta1)*1e3_dp, " ms"
    write(*,'(a,i4,a)') "Number of members in set: 1 + ", set%max_mem

    call getQ2min(0,Q2minPDF)
    call getQ2max(0,Q2maxPDF)
    set%Qmin = sqrt(Q2minPDF)
    set%Qmax = sqrt(Q2maxPDF)
    call getxmin(0,set%xmin)
    call getxmax(0,set%xmax)
    call getorderas(set%orderPDF) ! NB: LHAPDF returns 0 for 1-loop running, 1 for 2-loop etc.
    set%nloop = 1 + set%orderPDF
    call getthreshold(4,set%mc)
    call getthreshold(5,set%mb)
    call getthreshold(6,set%mt)
    call getnf(set%nfmax)
    if (set%nfmax .lt. 6) set%mt = 2.0d0*set%Qmax ! If no top in PDF set threshold beyond table max

    write(*,'(a,i4,a)') bold//"LHAPDF set: "//trim(name)//", member ",imem," loaded successfully"//reset

    ! Now let us define some hoppet specific parameters. These are
    ! typical values, and should guarantee similar accuracy as can be
    ! expected from LHAPDF
    ymax = real(ceiling(-log(set%xmin)), kind=dp) ! To get a nice value of ymax that can contain the full LHAPDF grid
    dy = 0.05_dp
    dlnlnQ = dy/4.0_dp
    if(ymax > 15.0) dlnlnQ = dy/8.0_dp ! for large ymax we need a finer grid in Q
    order = -6 ! Default
    yorder = 2 ! Quadratic interpolation in y
    lnlnQorder = 2 ! Quadratic interpolation in lnlnQ

    call hoppetSetPoleMassVFN(set%mc,set%mb,set%mt) ! set the pole masses
    call hoppetSetYLnlnQInterpOrders(yorder, lnlnQorder) ! Set the interpolation orders
    call cpu_time(ta1)
    call hoppetStartExtended(ymax, dy, set%Qmin, set%Qmax, dlnlnQ, set%nloop, order, factscheme_MSbar) ! Start hoppet
    call cpu_time(ta2)
    write(*,tfmt) "Time to start HOPPET: ", (ta2-ta1)*1e3_dp, " ms"
    write(*,'(a)') "Hoppet started with:"
    write(*,*) " ymax:       ", ymax
    write(*,*) " dy:         ", dy
    write(*,*) " Qmin:       ", set%Qmin
    write(*,*) " Qmax:       ", set%Qmax
    write(*,*) " dlnlnQ:     ", dlnlnQ
    write(*,*) " nloop:      ", set%nloop
    write(*,*) " order:      ", order
    write(*,*) " yorder:     ", yorder
    write(*,*) " lnlnQorder: ", lnlnQorder

    ! Now we fill the hoppet grid using the LHAPDF grid directly,
    ! rather than evolving ourselves
    !Q0 = set%Qmin
    call hoppetSetCoupling(alphasPDF(mz), mz, set%nloop)
    call cpu_time(ta1)
    call hoppetAssign(EvolvePDF)
    call cpu_time(ta2)
    write(*,tfmt) "Time to fill HOPPET grid from LHAPDF: ", (ta2-ta1)*1e3_dp, " ms"
    write(*,*) ! a blank line for clarity

    set%tables(0:0) => tables(0:0)
    set%owns_data = .false.

  end subroutine LoadLHAPDF

  !---------------------------------------------------
  subroutine LoadLHAPDFSet(name)    
    character(len=*),  intent(in) :: name
    !---------------------------------------
    type(hoppet_lhapdf_set), pointer :: set
    real(dp) :: ta1, ta2
    integer :: imem

    set => lhapdf_set ! shorthand

    ! load member 0; this also starts up Hoppet and sets the global tables(0) object
    call LoadLHAPDF(name, 0)

    allocate(set%tables(0:set%max_mem))
    call AllocPdfTable(set%tables, tables(0)) ! allocate the tables based on tables(0) structure    
    set%owns_data = .true.


    write(*,'(a)') bold//"Loading "//to_string(set%max_mem+1)//&
                   " members into hoppet from LHAPDF set: "//trim(set%name)//reset
    call cpu_time(ta1)
    do imem = 0, set%max_mem
      call initPDF(imem)
      call FillPdfTable_LHAPDF(set%tables(imem), EvolvePDF)
    end do
    call cpu_time(ta2)
    write(*,tfmt) "Time to load and fill all members into HOPPET: ", (ta2-ta1)*1e3_dp, " ms"
    write(*,*) !
  end subroutine LoadLHAPDFSet
end module hoppet_lhapdf