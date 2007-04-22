#!/usr/bin/perl -w

$basicopts="../small_fast_tab -nrep 1 -nxQ 5000  -Qmax 1e4 -ymax 11.5 -output";

$highprec="-dy 0.025 -dlnlnQ 0.005 -du 0.005 -olnlnQ 4 -order 6 -nopreev";


$highQ = "-dlnlnQ 0.005 -du 0.005 -olnlnQ 4";
$highY = "-dy 0.025 -order 6 -nopreev";
$highprec = "$highQ $highY";


# @dyvals=(0.04,0.05,0.065,0.08,0.1,0.15,0.2,0.25,0.3,0.4);
# 
# foreach $dy (@dyvals) {
#   $nameQ="dlnlnQ0.005-du0.005-olnlnQ4";
#   $nameY="dy$dy-order5-hires15-nopre";
#   $name="$nameY-$nameQ.res";
#   $command="$basicopts $highQ -dy $dy -order 5 -hires 15 -nopreev > res/$name";
#   print $command."\n";
#   system("$command")
# }

# @dlvals=(0.008,0.01,0.015,0.025,0.04,0.06,0.08,0.10,0.15);
# 
# foreach $dl (@dlvals) {
#   $nameQ="dlnlnQ$dl-du0.005-olnlnQ4";
#   #$nameY="dy0.025-order6-nopre";
#   $nameY="dy0.025-order6";
#   $highY = "-dy 0.025 -order 6";
#   $name="$nameY-$nameQ.res";
#   $command="$basicopts -dlnlnQ $dl -du 0.005 -olnlnQ 4 $highY > res/$name";
#   print $command."\n";
#   system("$command")
# }


# for time v. acc
@dyvals=(0.04,0.05,0.065,0.08,0.1,0.15,0.2,0.25,0.3,0.4);

@olnlnQ=(3,4);
@ord=(6,-6,5,-5);

foreach $olnlnQ (@olnlnQ) {
foreach $ord (@ord) {
foreach $dy (@dyvals) {
  $dl=$dy/4; # seems like reasonable choicie
  $argsQ="-dlnlnQ $dl -du 0.4 -olnlnQ $olnlnQ";
  $argsY="-dy $dy -order $ord -hires 15";
  #$argsY="-dy $dy -order $ord -hires 15 -nopreev";
  
  ($name="$argsY $argsQ.res") =~ s/(^-| )//g;
  $command="$basicopts $argsY $argsQ > res-pair/$name";
  print $command."\n";
  system("$command");
  $nrep=int(2e6*$dl**3);
  $command =~ s/-nxQ 5000//;
  $command =~ s/-nrep 1/-nrep $nrep/;
  $command =~ s/\.res/-speed.res/;
  print $command."\n";
  system("$command");
}
}
}
