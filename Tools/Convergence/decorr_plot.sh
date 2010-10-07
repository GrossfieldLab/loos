#!/bin/bash -x
#
# Basic plotting of decorr_time results (using defaults, mostly)
#
# Usage- decorr_plot.sh datafile title scale
#


#  This file is part of LOOS.
#
#  LOOS (Lightweight Object-Oriented Structure library)
#  Copyright (c) 2010, Tod D. Romo
#  Department of Biochemistry and Biophysics
#  School of Medicine & Dentistry, University of Rochester
#
#  This package (LOOS) is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation under version 3 of the License.
#
#  This package is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.


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
