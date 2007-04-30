#!/usr/bin/perl -w

# user option -go means we acutally run...
$go=($#ARGV >= 0 && $ARGV[0] =~ /-go/);

$basicopts="../prec_and_timing -nrep 1 -nxQ 5000  -Qmax 1e4 -ymax 11.5 -outputgrid";


$highQ = "-dlnlnQ 0.005 -du 0.005 -olnlnQ 4";
$highY = "-dy 0.025 -order 6 -4grids -nopreev";
#$highY15 = "-dy 0.025 -order 6 -hires15 -nopreev";
#$highprec = "$highY $highQ";


# # reference run
#runit("$highY $highQ -eps 1e-8");
#runit("$highY $highQ -asdt 0.1");


# # run for getting dy dependence
# @dyvals=(0.04,0.05,0.065,0.08,0.1,0.15,0.2,0.25,0.3,0.4);
# @opts=("-4grids -nopreev","-hires 15 -nopreev","-nopreev","-4grids","-hires 15","");
# foreach $opt (@opts) {
# foreach $dy (@dyvals) {
#   #$argsY="-dy $dy -order -6 -nopreev";
#   $argsY="-dy $dy -order -6 $opt";
#   runit("$argsY $highQ");
# }
# }

# run for getting dQ dependence

@dlvals=(0.008,0.01,0.015,0.025,0.04,0.06,0.08,0.10,0.15);
@opts=("-olnlnQ 3","-olnlnQ 4");
foreach $opt (@opts) {
foreach $dl (@dlvals) {
  $argsQ="-dlnlnQ $dl -du 0.005 $opt";
  runit("$highY $argsQ");
}
}

# # for time v. acc
# @dyvals=(0.04,0.05,0.065,0.08,0.1,0.15,0.2,0.25,0.3,0.4);
# 
# #@olnlnQ=(3,4);
# @olnlnQ=(4);
# @ord=(6,-6,-5);
# 
# foreach $olnlnQ (@olnlnQ) {
# foreach $ord (@ord) {
# foreach $dy (@dyvals) {
#   $dl=$dy/4; # seems like reasonable choicie
#   $argsQ="-dlnlnQ $dl -du 0.4 -olnlnQ $olnlnQ";
#   #$argsY="-dy $dy -order $ord";
#   #$argsY="-dy $dy -order $ord -hires 15 -nopreev";
#   $argsY="-dy $dy -order $ord -4grids -nopreev";
# 
#   $nrep=int(2e6*$dl**3);
#   runit("$argsY $argsQ", $nrep);
# }
# }
# }
# 

# build a name from a bunch of args
sub getname {
  ($args) = @_;
  ($name="$args.res") =~ s/(^-| )//g;
  return $name;
}

# run the program with the supplied args (basicopts assumed fixed)
# runit (args, nrep)
sub runit {
  @args = @_;
  $args = $args[0];
  $name = getname($args);
  $command="$basicopts $args > res-grid/$name";
  print $command."\n";
  if ($go) {system("$command")}
  # option test of speed
  if ($#args>0) {
    $nrep = $args[1];
    $command =~ s/-nxQ 5000//;
    $command =~ s/-nrep 1/-nrep $nrep/;
    $command =~ s/\.res/-speed.res/;
    print $command."\n";
    if ($go) {system("$command");}
  }
}

