#!/usr/bin/env bash

for i in `seq 7 18`
do
    flutile trim ha1 --subtype=H${i} --conversion=dna2aa H${i}.fna | smof clean -t p -drux > .H${i}
    N=`smof wc -l .H${i}`
    for j in `seq 1 $N`
    do
        diff <(smof cut -f $j H${i}.ha1 | smof stat) <(smof cut -f $j .H${i} | smof stat) 2> /dev/null
        if [ $? -ne 0 ]
        then

            echo " -- HA1 trimming failed on test H${i} entry ${j}"
            echo " -- Alignment of GenBank HA1 versus flutile:"
            cat <(smof cut -f $j .H${i}) <(smof cut -f $j H${i}.ha1) |
                mafft --quiet --clustalout /dev/stdin
            exit 1
        fi
    done
done
