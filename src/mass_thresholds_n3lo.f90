!======================================================================
! What follows is stuff for VFNS at order(as^3)
!
! As of 2024-01, the list of relevant articles is to be found at
! https://github.com/MSHTPDF/N3LO_additions
! 
!  J. McGowan, et. al.,    arXiv:2207.04739 -- MSHT parametrisations
!  H. Kawamura, et. al.,   arXiv:1205.5727  -- KPMV parametrisations (A3PShg & A3PShq, eqs 3.49-3.53, old moments)
!  I. Bierenbaum, et. al., arXiv:0904.3563  -- BBK: moments up to N=12-14, for all structures, defs in (6.1-6.3)
!  J. Ablinger, et. al.,   arXiv:1406.4654  -- A3NSqq_H, all moments & x-space (App.B, even/odd to be understood)
!  J. Ablinger, et. al.,   arXiv:1409.1135  -- A3PShq, all moments & x-space (Eq. 5.41)
!  J. Ablinger, et. al.,   arXiv:1402.0359  -- A3Sgq_H, x-space (Eq 6.45)
! 
!  J. Ablinger, et. al.,   arXiv:2211.05462 -- A3Sgg_H
!
!
! Pieces that are going to be needed for my purposes are those from that
! code which are free of logs of mu^2/m^2, since we will for the time being
! keep these two scales equal (and later on, if we do not then pieces
! with logs will be reconstructed appropriately?)
!
! Thus the pieces needed, in the notation of the second of the above 
! papers, are:
!
! For Delta (f_k + f_kbar):
!   A^{3,NS}_{qq,H}    A3NSqq_H
!   A^{3,PS}_{qq,H}    A3PSqq_H
!   A^{3,S}_{qg,H}     A3Sqg_H
!
! For Delta (f_H + f_Hbar)
!   A^{3,PS}_{Hq}      A3PShg
!   A^{3,S}_{Hg}       A3PShq
!
! For Delta (G)
!   A^{3,S}_{gq,H}     A3Sgq_H
!   A^{3,S}_{gg,H}     A3Sgg_H
!
!----------------------------------------------------------------------

module mass_thresholds_n3lo
  use types; use consts_dp; use convolution_communicator
  implicit none
  private


end module mass_thresholds_n3lo
