#!/bin/bash

#set -x

function run {
	
	k=$1 s=$2 D=$3 t=$4 I=$5 U=$6 o=$7 T=$8 w=$9 e=${10} x=${11}
	
	if [ $U == 6.7 ] ; then
		UU=0670
	else
		UU=1788
	fi
	
	filename=LSout_s${s}_k${k}_U${UU}_I${I}_o${o}_b
    args="-k${k}0 -s${s} -I${I}000 -U${U} -o${o} -T${T} -w${w} -e${e} -f${filename} -x${x}"
    
    echo "executing: ./LSfire+ "${args}
	./LSfire+ ${args}
	echo ${filename} done
	./move ${filename}
}
x=$1  # current time step from ForeFire
D=25  # diffusivity (gaussian)
t=1   # heating-before-burning model parameter
U=3 # wind velocity [m/s]
I=20  # fire intensity [MW/m]
w=1   # constant wind
e=10  # save each 10 iterations
s=1   # basic grid

##################### lsm | without/with obstacles #####################
#x=$1 k=0 o=0 T=300 ; run $k $s $D $t $I $U $o $T $w $e $f $x #  *** OK ***  tmax: 181.1
#k=0 o=2 T=300 ; run $k $s $D $t $I $U $o $T $w $e 

################## gaussian | without/with obstacles ##################
#k=1 o=0 T=300 ; run $k $s $D $t $I $U $o $T $w $e #  *** OK ***  tmax: 138.5
x=$1 k=5 o=1 T=300 ; run $k $s $D $t $I $U $o $T $w $e $f $x #  *** OK ***  tmax: Tmax

############# gaussian*lognorm | without/withn obstacles ################
#k=3 o=0 T=300 ; run $k $s $D $t $I $U $o $T $w $e # *** OK *** tmax: 101.5
#x=$1 k=6 o=1 T=300 ; run $k $s $D $t $I $U $o $T $w $e $f $x # *** OK *** tmax: 151.0

