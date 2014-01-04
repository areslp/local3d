#!/bin/sh
#for file in `ls *.obj| grep -v noise`
#do
	file=box.obj
    echo $file
    ./addnoise.exe $file 0.1
    ./addnoise.exe $file 0.2
    ./addnoise.exe $file 0.3
    ./addnoise.exe $file 0.4
    ./addnoise.exe $file 0.5
#done
