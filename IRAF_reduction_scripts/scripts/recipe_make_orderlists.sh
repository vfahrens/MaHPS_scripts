#!/bin/bash

fname="*_ods_fred.fits"
for i in {1..84..1}
do
	let k=i+64
  rm fxcor_ord$k.lis
  for j in `echo $fname`
    do
	    echo $j[$i] >> fxcor_ord$k.lis
    done
done
