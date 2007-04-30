#!/usr/bin/perl -w

$command=$0." ".join(" ",@ARGV);
print "# $command\n";

$runfile=$ARGV[0];

#$reffile="res/dy0.025-order6-nopre-dlnlnQ0.005-du0.005-olnlnQ4.res";
$reffile="res-grid/reference.res";
#$reffile=$ARGV[1];

$refls=`ls -l $reffile`;
print "# Reffile: $refls";


($runfileA = $runfile ) =~ s/(.*)XX.*/$1/;
($runfileB = $runfile ) =~ s/.*XX(.*)/$1/;

($dir,$runfileA) = split("/",$runfileA);

open(LS,"ls $dir/|") || die "could not get listing";
while ($file=<LS>) {
  if ($file =~ /^$runfileA([0-9.]+)$runfileB$/) {
    chomp($file);
    $val = $1;
    print "$val ";
    $command="./compare2files_v2 $dir/$file $reffile -summary -protect";
    #print "$command\n";
    $res=`$command`;
    chomp $res;
    print "$res  # $command\n";
  }
}


