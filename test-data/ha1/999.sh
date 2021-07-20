#!/usr/bin/env bash

for i in `seq 1 18`
do
    sed "s/__HA_SUBTYPE__/H${i}/" 001.rq > 001-H${i}.rq
    octofludb query 001-H${i}.rq | sort -R | head -50 > 100-H${i}.txt

    cat 100-H${i}.txt |
        awk 'BEGIN {FS="\t"; OFS="\n"} {print ">" $1 "|" $2 "|" $3 , $4}' | smof clean -t n -dux > H${i}.fna
    cat 100-H${i}.txt |
        awk 'BEGIN {FS="\t"; OFS="\n"} {print ">" $1 "|" $2 "|" $3 , $5}' | smof clean -t p -drux > H${i}.faa
    cat 100-H${i}.txt |
        awk 'BEGIN {FS="\t"; OFS="\n"} {print ">" $1 "|" $2 "|" $3 , $6}' | smof clean -t p -drux > H${i}.ha1
done
