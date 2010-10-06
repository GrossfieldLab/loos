#!/bin/bash -x
#
# (c) 2010 Tod D. Romo, Grossfield Lab, URMC
#
# Basic plotting of decorr_time results (using defaults, mostly)
#
# Usage- decorr_plot.sh datafile title scale
#

DATA="$1"
BASE=`echo $DATA | cut -d. -f1`
TITLE=${2:-Decorrelation Time Analysis}
SCALE=${3:-1}


RANGE=`perl -ane 'BEGIN{$min=1e100;$max=0;}{next if(/^#/);if($F[0]>$max){$max=$F[0];}if($F[0]<$min){$min=$F[0];}}END{print"[$min:$max]";}' $DATA`
echo "$RANGE"

gnuplot <<EOF
set out "$BASE.ps"
set term post enhanced solid color landscape
set log x
set title "$TITLE"
set xlabel '{/Symbol t}_d'
set ylabel '{/Symbol s}^2'
set xrange $RANGE

f(x) = 1

plot "$DATA" u (\$1/$SCALE):2:3 w errorl ti 'N=2',\
'' u (\$1/$SCALE):4:5 w errorl ti 'N=4',\
'' u (\$1/$SCALE):6:7 w errorl ti 'N=10',\
'' u (\$1/$SCALE):(f(\$1/$SCALE)) w l lw 4 not

EOF


ps2pdf $BASE.ps $BASE.pdf ; rm $BASE.ps
