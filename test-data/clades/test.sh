#!/usr/bin/env bash

set -u
set -e
set -v

# these should all yield the same final amino acid sequences (gaps are removed
# since `smof translate` does not preserve gap spacing)
flutile trim ha1 --subtype=H1 --conversion=dna2aa  1A.1-like.fna | smof clean -xu > .a
flutile trim ha1 --subtype=H1 --conversion=dna2dna 1A.1-like.fna | smof translate | smof clean -xu > .b
flutile trim ha1 --subtype=H1 --conversion=aa2aa   1A.1-like.faa | smof clean -xu > .c
diff .a .b
diff .b .c


# same for h3
flutile trim ha1 --subtype=H3 --conversion=dna2aa  3.2010.2.fna | smof clean -xu > .a
flutile trim ha1 --subtype=H3 --conversion=dna2dna 3.2010.2.fna | smof translate | smof clean -xu > .b
flutile trim ha1 --subtype=H3 --conversion=aa2aa   3.2010.2.faa | smof clean -xu > .c
diff .a .b
diff .b .c

rm -f .a .b .c
