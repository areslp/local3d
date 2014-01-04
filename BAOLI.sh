#!/bin/sh
#file_list='bunny double-torus1 fandisk2 mechpart'
#noise_level='0.1 0.2 0.3 0.4 0.5'
file_list='box'
noise_level='0.2 0.5'
cd baoli
for file in $file_list
do
    echo "==========================================================================="
    echo $file "BAOLI noise level 0" 
    #baoli 
    ./normal_estimation.exe ../paper_data/$file.xyz 30 $file.xyzn
    ./xyzn2apts $file.xyzn
    time ./RobustNormalEstimationCMD.exe $file.apts out.apts -n -o 0.7
    ./apts2xyzn out.apts
    matlab -nodesktop -nosplash -r -wait "rearrange_points('out.xyzn','../paper_data/ori_data/$file.xyzn');quit;"
    ./consistent_normal.exe out.xyzn --ori "../paper_data/ori_data/$file.xyzn"

    for nl in $noise_level
    do
        echo $file "BAOLI noise level " $nl
        ./normal_estimation.exe "../paper_data/noise_"$nl"_"$file".xyz" 30 $file.xyzn
        ./xyzn2apts $file.xyzn
        time ./RobustNormalEstimationCMD.exe $file.apts out.apts -n -o 0.7
        ./apts2xyzn out.apts
        matlab -nodesktop -nosplash -r -wait "rearrange_points('out.xyzn','../paper_data/ori_data/$file.xyzn');quit;"
        ./consistent_normal.exe out.xyzn --ori "../paper_data/ori_data/$file.xyzn"
    done

done
