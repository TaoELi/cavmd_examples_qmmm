#!/bin/bash

split -l 30580 simu_1.xc.xyz --verbose tmp-


for f in tmp*
do
    echo "doing $f"
    cp $f $f.xyz
    rm $f
    # now each contains 16 different local geometries
    sed -i '30579,30580d' $f.xyz
    sed -i '1,2d' $f.xyz
    # new file name
    fnew="$f-"
    split -l 1911 $f.xyz --verbose $fnew

    rm $f.xyz
    
    for f2 in $fnew*
    do
	mv $f2 $f2.xyz
	sed -i "1 i\1911\n" $f2.xyz
    done
done

rm 1.xyz 2.xyz 3.xyz 4.xyz 5.xyz 6.xyz 7.xyz 8.xyz 9.xyz 10.xyz 11.xyz 12.xyz 13.xyz 14.xyz 15.xyz 16.xyz
for f in tmp-*-aa.xyz
do
    sed -i "4s/C/F/" $f
    sed -i "5s/C/F/" $f
    cat $f >> 1.xyz
done
for f in tmp-*-ab.xyz
do
    sed -i "4s/C/F/" $f
    sed -i "5s/C/F/" $f
    cat $f >> 2.xyz
done
for f in tmp-*-ac.xyz
do
    sed -i "4s/C/F/" $f
    sed -i "5s/C/F/" $f
    cat $f >> 3.xyz
done
for f in tmp-*-ad.xyz
do
    sed -i "4s/C/F/" $f
    sed -i "5s/C/F/" $f
    cat $f >> 4.xyz
done
for f in tmp-*-ae.xyz
do
    sed -i "4s/C/F/" $f
    sed -i "5s/C/F/" $f
    cat $f >> 5.xyz
done
for f in tmp-*-af.xyz
do
    sed -i "4s/C/F/" $f
    sed -i "5s/C/F/" $f
    cat $f >> 6.xyz
done
for f in tmp-*-ag.xyz
do
    sed -i "4s/C/F/" $f
    sed -i "5s/C/F/" $f
    cat $f >> 7.xyz
done
for f in tmp-*-ah.xyz
do
    sed -i "4s/C/F/" $f
    sed -i "5s/C/F/" $f
    cat $f >> 8.xyz
done
for f in tmp-*-ai.xyz
do
    sed -i "4s/C/F/" $f
    sed -i "5s/C/F/" $f
    cat $f >> 9.xyz
done
for f in tmp-*-aj.xyz
do
    sed -i "4s/C/F/" $f
    sed -i "5s/C/F/" $f
    cat $f >> 10.xyz
done
for f in tmp-*-ak.xyz
do
    sed -i "4s/C/F/" $f
    sed -i "5s/C/F/" $f
    cat $f >> 11.xyz
done
for f in tmp-*-al.xyz
do
    sed -i "4s/C/F/" $f
    sed -i "5s/C/F/" $f
    cat $f >> 12.xyz
done
for f in tmp-*-am.xyz
do
    sed -i "4s/C/F/" $f
    sed -i "5s/C/F/" $f
    cat $f >> 13.xyz
done
for f in tmp-*-an.xyz
do
    sed -i "4s/C/F/" $f
    sed -i "5s/C/F/" $f
    cat $f >> 14.xyz
done
for f in tmp-*-ao.xyz
do
    sed -i "4s/C/F/" $f
    sed -i "5s/C/F/" $f
    cat $f >> 15.xyz
done
for f in tmp-*-ap.xyz
do
    sed -i "4s/C/F/" $f
    sed -i "5s/C/F/" $f
    cat $f >> 16.xyz
done

rm tmp*.xyz
