#!/bin/bash


#set -x

function run {
	
	k=$1 s=$2 D=$3 t=$4 I=$5 U=$6 o=$7 T=$8 w=$9 e=${10}
	
	if [ $U == 6.7 ] ; then
		UU=0670
	else
		UU=1788
	fi
	
	filename=LSout_s${s}_k${k}_U${UU}_I${I}_o${o}_b
    args="-k${k}0 -s${s} -I${I}000 -U${U} -o${o} -T${T} -w${w} -e${e} -f${filename}"
    
    echo "executing: ./LSfire+ "${args}
	./LSfire+ ${args}
	echo ${filename} done
	./move ${filename}
}

D=25  # diffusivity (gaussian)
t=1   # heating-before-burning model parameter
U=6.7 # wind velocity [m/s]
I=10  # fire intensity [MW/m]
w=1   # constant wind
e=10  # save each 10 iterations
s=1   # basic grid

##################### lsm | without/with obstacles #####################
k=0 o=0 T=200 ; run $k $s $D $t $I $U $o $T $w $e #  *** OK ***  tmax: 181.1
k=0 o=1 T=300 ; run $k $s $D $t $I $U $o $T $w $e 

################## gaussian | without/with obstacles ##################
k=1 o=0 T=200 ; run $k $s $D $t $I $U $o $T $w $e #  *** OK ***  tmax: 138.5
k=1 o=1 T=200 ; run $k $s $D $t $I $U $o $T $w $e #  *** OK ***  tmax: Tmax

############# gaussian*lognorm | without/withn obstacles ################
k=3 o=0 T=200 ; run $k $s $D $t $I $U $o $T $w $e # *** OK *** tmax: 101.5
k=3 o=1 T=200 ; run $k $s $D $t $I $U $o $T $w $e # *** OK *** tmax: 151.0

