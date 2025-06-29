
# our rescaling function

zeta_of_y(y)=-(y + 9*(1-exp(-y)))
zeta_of_x(x)=-(log(1.0/x) + 9*(1-x))

unset log x
set xrange [zeta_of_x(1e-5):0]
set xlabel 'x'

set xtics mirror \
         ("10^{-5}" zeta_of_x(1e-5),\
          "10^{-4}" zeta_of_x(1e-4),\
          "10^{-3}" zeta_of_x(1e-3),\
          "10^{-2}" zeta_of_x(1e-2),\
          "0.1" zeta_of_x(0.1),\
          ".2"  zeta_of_x(0.2) 1,\
          ".3"  zeta_of_x(0.3),\
          ".4"  zeta_of_x(0.4) 1,\
          ".5"  zeta_of_x(0.5),\
          ".6"  zeta_of_x(0.6) 1,\
          ".7"  zeta_of_x(0.7),\
          ".8"  zeta_of_x(0.8) 1,\
          ".9"  zeta_of_x(0.9),\
          "1"   zeta_of_x(1.0)\
          )



set xtics mirror  add \
          (" " zeta_of_x(2e-5) 1,\
           " " zeta_of_x(5e-5) 1,\
           " " zeta_of_x(2e-4) 1,\
           " " zeta_of_x(5e-4) 1,\
           " " zeta_of_x(2e-3) 1,\
           " " zeta_of_x(5e-3) 1,\
           " " zeta_of_x(2e-2) 1,\
           " " zeta_of_x(5e-2) 1,\
           " " zeta_of_x(2e-1) 1,\
           " " zeta_of_x(5e-1) 1\
          )

