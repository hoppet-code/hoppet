#!/usr/bin/perl -w

$command=$0." ".join(" ",@ARGV);
print "# $command\n";

$runfile=$ARGV[0];

#$reffile="res/dy0.025-order6-nopre-dlnlnQ0.005-du0.005-olnlnQ4.res";
$reffile="res-grid/reference.res";
#$reffile=$ARGV[1];

$refls=`ls -l $reffile`;
print "# Reffile: $refls";

($runfileA = $runfile ) =~ s/(.*)XX.*XX.*/$1/;
($runfileB = $runfile ) =~ s/.*XX(.*)XX.*/$1/;
($runfileC = $runfile ) =~ s/.*XX.*XX(.*)/$1/;
#$runfileA =~ s/res-pair\///;
#$reffile =~ s/res\///;

($dir,$runfileA) = split("/",$runfileA);

# remove the speed prefix if it was there (we'll put it back in later)
$runfileC=~ s/-speed//;

open(LS,"ls $dir/|") || die "could not get listing";
print "# format: dy dlnlnQ guds0.7 gudsc0.7 guds0.9 gudsc0.9 time-init time-preev time-ev # command\n";
while ($file=<LS>) {
  if ($file =~ /^$runfileA([0-9.]+)$runfileB([0-9.]+)$runfileC$/) {
    chomp($file);
    $dyval = $1;
    $dlval = $2;
    print "$dyval $dlval";
    $command="./compare2files_v2 $dir/$file $reffile -protect -summary";
    #print "$command\n";
    $res=`$command`;
    chomp $res;
    print "$res ";
    @last=split(/\s+/,`tail -1 $dir/$file`);
    print "$last[6] $last[7] ";
    $file =~ s/\.res/-speed.res/;
    @last=split(/\s+/,`tail -1 $dir/$file`);
    print "$last[8] ";

    print "  # $command\n";
  }
}


