#!/usr/bin/perl -w

$runfile=$ARGV[0];

$reffile="res/dy0.025-order6-nopre-dlnlnQ0.005-du0.005-olnlnQ4.res";
#$reffile=$ARGV[1];


($runfileA = $runfile ) =~ s/(.*)XX.*/$1/;
($runfileB = $runfile ) =~ s/.*XX(.*)/$1/;
$runfileA =~ s/res\///;
#$reffile =~ s/res\///;

open(LS,"ls res/|") || die "could not get listing";
while ($file=<LS>) {
  if ($file =~ /^$runfileA([0-9.]+)$runfileB$/) {
    chomp($file);
    $val = $1;
    print "$val ";
    $command="./compare2files res/$file $reffile";
    #print "$command\n";
    $res=`$command`;
    chomp $res;
    print "$res  # $command\n";
  }
}


