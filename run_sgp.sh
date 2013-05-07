#!/bin/sh
if [ $# -lt 1 ] 
then
    echo "./run_sgp.sh fname"
	exit 1
fi
fname=$1
k=200
if [ $# == 2 ]
then
	k=$2
fi
echo ${fname}
./normal_sgp.exe data/${fname}.xyz ${k}
./consistent_normal.exe sgp_out.xyzn --ori data/standard_normal.xyzn
