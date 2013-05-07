#!/bin/sh
if [ $# -lt 1 ] 
then
    echo "./run_sgp.sh fname"
	exit 1
fi
fname=$1
k=500
if [ $# == 2 ]
then
	k=$2
fi
echo ${fname}
echo ${k}
cp data/${fname}.xyzn baoli/
cd baoli
./xyzn2apts ${fname}.xyzn
./RobustNormalEstimationCMD.exe ${fname}.apts out.apts -n -o 0.7 -l ${k}
./apts2xyzn out.apts
cd ..
echo "plz run rearrange_points in matlab!"
read -n 1 -p "Press any key to continue..."
#matlab -nodesktop -nosplash -r "rearrange_points;quit;"
cd baoli
./consistent_normal.exe out.xyzn --ori ../data/standard_normal.xyzn
