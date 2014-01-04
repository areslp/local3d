#!/bin/sh
#file_list='bunny double-torus1 fandisk2 mechpart'
#noise_level='0.1 0.2 0.3 0.4 0.5'
file_list='box'
noise_level='0.2 0.5'
for file in $file_list
do
    echo "==========================================================================="
    echo $file "RRHT noise level 0" 
    # rrht
    time ./normal_sgp.exe paper_data/$file.xyz
	#matlab -nodesktop -nosplash -r -wait "rearrange_points('out_sgp.xyzn','../paper_data/ori_data/$file.xyzn');quit;"
    ./consistent_normal.exe sgp_out.xyzn --ori paper_data/ori_data/$file.xyzn

    for nl in $noise_level
    do
        echo $file "RRHT noise level " $nl
        time ./normal_sgp.exe "paper_data/noise_"$nl"_"$file".xyz"
		#matlab -nodesktop -nosplash -r -wait "rearrange_points('out_sgp.xyzn','../paper_data/ori_data/$file.xyzn');quit;"
        ./consistent_normal.exe sgp_out.xyzn --ori paper_data/ori_data/$file.xyzn
    done

    # echo "==========================================================================="
    # echo $file "BaoLi"
    # baoli
done
