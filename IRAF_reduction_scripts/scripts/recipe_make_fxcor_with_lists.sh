#!/bin/bash
name1="fxcor_ord"
name2=".lis"
for i in {065..148..1}
do
  echo $i
  for j in `ls $name1$i$name2`
  do
    echo "fxcor @"$j "20190903_0104_FOC1903_SCC2_ods_ord"$i"_A_lin_IRAF.fits" >> fxcor_with_lists.cl
  done
done