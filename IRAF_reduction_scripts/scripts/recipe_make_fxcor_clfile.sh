for i in `ls *.fits`;
    do echo "fxcor" $i $reffile${i:34} "width=20 maxwidth=20 interac=no autowri=yes output="${i:0:13}"_fxcor_result" >> fxcor_clfile.txt;
done
