if test x$donotrun != xtrue; then
    tmp=`mktemp -d`
    TIMESTEP="2e-3 1e-3 5e-4"

    for timestep in $TIMESTEP; do
        rm -rf k.dat
        if sed -e "s/FILTER/2/g" -e "s/TIMESTEP/$timestep/g" < $1 | gerris2D -; then :
            k0=$(awk 'NR == 1 {print $5}' k.dat)
            awk '{print $3, sqrt(($5/'$k0'-1)*($5/'$k0'-1))}'  k.dat > k-$timestep
	else
	    exit 1
	fi
    done
    rm -rf $tmp
fi

if cat <<EOF | gnuplot ; then :
  set term postscript eps enhanced color "Helvetica" 14
  set output "k-skew.eps"
  set ylabel '(k-k_0)/k_0'
  set xlabel 't'
  set yrange[0.:0.5]
  p "k-2e-3" u 1:2 title 'dt=2d-3' w l,  "k-1e-3" u 1:2 title '1e-3' w l, "k-5e-4" u 1:2 title '5e-4' w l 
EOF
else
    exit 1
fi

if cat <<EOF | python ; then :
from check import *
from sys import *
for name in ['k-2e-3','k-1e-3','k-5e-4']:
  if Curve(name,2,2).max() > 0.05:
    exit(1)
EOF
else
   exit 1
fi
