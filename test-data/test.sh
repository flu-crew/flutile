#!/usr/bin/env bash

set -e
set -u

diff <(flutile aadiff --subtype=H1 h1-aadiff-ex.faa) h1-aadiff-ex.txt

flutile aadiff --caton82 --annotation-tables=rbs.txt --subtype=H1 --join-annotations 1B.2.1-aadiff-test.faa > .obs-aadiff
diff .obs-aadiff .exp-aadiff

flutile aadiff --wiley81 --subtype=H3 3.1990.1-aadiff.faa > .obs-h3-aadiff
diff .obs-h3-aadiff .exp-h3-aadiff

flutile annotate --subtype=H1 --caton82 --annotation-tables=rbs.txt --join-annotations  1B.2.1-aadiff-test.faa > .obs-annotate
diff .obs-annotate .exp-annotate

flutile trim ha1 --subtype=H1 1B.2.1.fna > .obs.ha1

# ha1 can also be extracts as a motif relative to the start of the mature protein
diff .obs.ha1 <(flutile trim motif --subtype=H1 -m "1-327" --fasta 1B.2.1.fna)

# or the motif can be numbered relative to initial methionine
diff .obs.ha1 <(flutile trim motif --subtype=H1 -m "18-344" --keep-signal --fasta 1B.2.1.fna)

flutile trim motif --subtype=H3 -m "motif=145,155,156,158,159,189" clade4.fna > .obs-h3-motif
diff .obs-h3-motif .exp-h3-motif

diff <(smof subseq -b 36 59 wsn33.faa) \
     <(cat N1.fna wsn33.fna | flutile trim motif --subtype=N1 -m "stalk=36-59" --fasta | smof grep WSN | smof clean -x)

for seq in `seq 1 18`
do
    smof sample --seed=42 --number=10 H${seq}.faa |
        flutile aadiff --subtype=H${seq} > .obs-H${seq}-aadiff.txt
    diff .obs-H${seq}-aadiff.txt .H${seq}-aadiff.txt
done

for seq in `seq 1 11`
do
    smof sample --seed=42 --number=10 N${seq}.fna |
        smof translate -sf |
        flutile aadiff --subtype=N${seq} > .obs-N${seq}-aadiff.txt
    diff .obs-N${seq}-aadiff.txt .N${seq}-aadiff.txt
done

for seq in PB1 PB2 PA NP M NS
do
    smof sample --seed=42 --number=10 ${seq}.fna |
        smof translate -sf |
        flutile aadiff --subtype=${seq} > .obs-${seq}-aadiff.txt
    diff .obs-${seq}-aadiff.txt .${seq}-aadiff.txt
done

rm .obs-*
