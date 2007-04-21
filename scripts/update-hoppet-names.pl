#!/usr/bin/perl -w
#======================================================================
# script to convert files from old naming schemes to new ones.
#
# Usage:
#  update-hoppet-names.pl file.f90
#
# If there are any names to be modified, the original file will be copied
# to file.f90.bak and file.f90 will be modified to contain new names.
#
# Note: if you are using this for the first time on you own programs
# check that nothing too strange has happened before throwing away the
# backup... Note two very short names are renamed: sh (->dh) and gd
# (-> grid) they could conceivably conflict with something a user
# program.
#
# Author: Gavin Salam, 2006
# License: GNU Public license
#======================================================================

# names that are probably OK, but could do with removal
# of the _conv prefix
@removeconvprefix = (
"InitGridDef",
"AllocGridQuant",
"InitGridQuant",
"InitGridQuantSub",
"InitGridQuantLHAPDF",
"PrintGridQuant",
"EvalGridQuant",
"MomGridQuant",
"WgtGridQuant",
"AllocGridConv",
"InitGridConv",
"ValidateGD",
"GetDerivedProbes",
"SetDerivedConv",
"SetDerivedConv_nodealloc",
);

@removepdfgenprefix = (
"AllocPDF",
"InitPDF",
"InitPDFSub",
"InitPDF_LHAPDF",
"AllocInitPDF",
"AllocInitPDFSub",
);

# names that need to made more sensible (in convolution.f90)
$rename{qr/conv\\?_DelGridQuant/i} = "Delete";
$rename{qr/conv\\?_DelGridConv/i}  = "Delete";
$rename{qr/conv\\?_ZeroGridConv/i} = "SetToZero";
$rename{qr/conv\\?_MultGridConv/i} = "Multiply";
$rename{qr/conv\\?_AddGridConv/i}  = "AddWithCoeff";
$rename{qr/conv\\?_ConvGridConv/i} = "SetToConvolution";
$rename{qr/conv\\?_CommGridConv/i} = "SetToCommutator";
$rename{qr/conv\\?_Seteps/i} = "SetConvolutionEps";
$rename{qr/SetConvolutionEps/i} = "SetDefaultConvolutionEps";



# names that need to made more sensible (in conv_objects)
$rename{qr/cobj\\?_InitSplitPolLO/i}   = "InitSplitMatPolLO";
$rename{qr/cobj\\?_InitSplitPolNLO/i}  = "InitSplitMatPolNLO";
$rename{qr/cobj\\?_InitSplitLO/i}   = "InitSplitMatLO";
$rename{qr/cobj\\?_InitSplitNLO/i}  = "InitSplitMatNLO";
$rename{qr/cobj\\?_InitSplitNNLO/i} = "InitSplitMatNNLO";
$rename{qr/cobj\\?_InitMTMNNLO/i} = "InitMTMNNLO";


# pdf representations
$rename{qr/pdfr\\?_LabelRep/i}   = "LabelPdfAsRep";
$rename{qr/pdfr\\?_GetRep/i}     = "GetPdfRep";
$rename{qr/pdfr_HumanToEvln/i}   = "CopyHumanPdfToEvln";
$rename{qr/pdfr_EvlnToHuman/i}   = "CopyEvlnToHumanPdf";
$rename{qr/CopyEvlnToHumanPdf/i} = "CopyEvlnPdfToHuman";# correct above mistake!

# conv_objects -> dglap_objects
$rename{qr/cobj\\?_DefaultPrep/i} = "DefaultEvlnRep";
$rename{qr/conv\\?_objects/i}     = "dglap_objects";
$rename{qr/PMat/i}                = "split_mat";
$rename{qr/CMat/i}                = "coeff_mat";
$rename{qr/MassThresholdMat/i}    = "mass_threshold_mat";
$rename{qr/cobj_InitSplit/i}      = "InitSplitMat";
# things we'll have to sort out partially separately (automated
# name change will deal with calls, but not definition and interface)
$rename{qr/cobj_AddSplit/i}       = "AddWithCoeff";
$rename{qr/cobj_DelSplit/i}       = "Delete";
$rename{qr/cobj_AllocSplit/i}     = "AllocSplitMat";
$rename{qr/cobj_SetNfMTM/i}       = "SetNfMTM";
$rename{qr/cobj\\?_GetDerivedProbes/i} = "GetDerivedSplitMatProbes";
$rename{qr/cobj_SetDerivedSplit/i}= "SetDerivedSplitMat";


# sigma_holder -> dglap_holder
$rename{qr/holders/i}             = "dglap_holders";
$rename{qr/sigma_holder/i}        = "dglap_holder";
$rename{qr/holder_InitSigma/i}    = "InitDglapHolder";
$rename{qr/holder_SetNf/i}        = "SetDglapHolderNf";
$rename{qr/SetDglapHolderNf/i}    = "SetNfDglapHolder";
$rename{qr/sh/i}                  = "dh";
$rename{qr/gd/i}                  = "grid";
$rename{qr/dh\%P2/i}              = "dh%P_NNLO";
$rename{qr/dh\%P1/i}              = "dh%P_NLO";
$rename{qr/dh\%P/i}               = "dh%P_LO";


# alphas
$rename{qr/as_Del/i}    = "Delete";
$rename{qr/as_server/i} = "qcd_coupling";
$rename{qr/as_handle/i} = "running_coupling";
$rename{qr/ash/i}       = "coupling";
$rename{qr/as_Init/i}   = "InitRunningCoupling";
$rename{qr/as_Value/i}  = "Value";
$rename{qr/as_nfRange/i}     = "NfRange";
$rename{qr/as_nfAtQ/i}       = "NfAtQ";
$rename{qr/as_QRangeAtNf/i}  = "QRangeAtNf";

# evolution
$rename{qr/ev_evolve_varnf/i}      = "EvolvePDF";
$rename{qr/ev_evolve_varnf_gen/i}  = "EvolveGeneric";
$rename{qr/ev_DelEvOp/i}           = "Delete";
$rename{qr/ev\\?_Setdu/i}            = "SetDefaultEvolutionDu";
$rename{qr/ev\\?_Setdt/i}            = "SetDefaultEvolutionDt";

# pdf tabulation
$rename{qr/pdf_tabulate_new/i} = "pdf_tabulate";
$rename{qr/pdftab_DelTab/i} = "Delete";
$rename{qr/pdftab/i} = "pdf_table";
$rename{qr/pdftab_AllocTab/i} = "AllocPdfTable";
$rename{qr/pdftab_AssocNfInfo/i} = "AddNfInfoToPdfTable";

$rename{qr/pdftab_InitTabSub/i    }  = "FillPdfTable_f90sub";
$rename{qr/pdftab_InitTab_LHAPDF/i}  = "FillPdfTable_LHAPDF";
$rename{qr/pdftab_InitTabEvolve/i }  = "EvolvePdfTable";
$rename{qr/pdftab_TabEvolveGen/i }   = "EvolvePdfTableGen";
$rename{qr/pdftab_PreEvolve/i     }  = "PreEvolvePdfTable";
$rename{qr/pdftab_ValTab_yQ/i     }  = "EvalPdfTable_yQ";
$rename{qr/pdftab_ValTab_xQ/i     }  = "EvalPdfTable_xQ";

# outside access
$rename{qr/pdfevln_all_public/i} = "hoppet_v1";

$result = '';

if ($#ARGV >= 0) {
  open ($filehandle, "<$ARGV[0]") || die "Could not open $ARGV[0] for reading";
} else {
  $filehandle = STDIN;
}

while ($line=<$filehandle>) {
  $result .= $line;
}
close($filehandle);

$orig_result = $result;

foreach $name (@removeconvprefix) {
  $prefixname = qr/conv\\?_$name/;
  $result =~ s/\b$prefixname\b/$name/gi;
}
foreach $name (@removepdfgenprefix) {
  $prefixname = qr/pdfgen\\?_$name/;
  $result =~ s/\b$prefixname\b/$name/gi;
}

foreach $name (keys %rename) {
  $newname = $rename{$name};
  $result =~ s/\b$name\b/$newname/gi;
}

# some more delicate cases
$result =~ s/\bid_([a-z]+)\b/iflv_$1/gi;


if ($#ARGV >= 0) {
  $name = $ARGV[0];
  if ($result eq $orig_result) {
    print STDERR "No modifications need be made to $name\n";
  } else {
    print STDOUT "Renaming $name to $name.bak\n";
    rename ($name, "$name.bak");
    open (OUT, ">$name") || die "Could not open $name";
    print STDOUT "Writing new version to $name\n";
    print OUT $result;
  }
} else {
  print $result;
}
