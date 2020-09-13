#!/usr/bin/env bash

set -u
set -e

# these should all yield the same final amino acid sequences (gaps are removed
# since `smof translate` does not preserve gap spacing)
flutile trim h1-ha1 --conversion dna2aa  1A.1-like.fna 2> /dev/null | smof clean -xu > .a
flutile trim h1-ha1 --conversion dna2dna 1A.1-like.fna 2> /dev/null | smof translate | smof clean -xu > .b
flutile trim h1-ha1 --conversion aa2aa   1A.1-like.faa 2> /dev/null | smof clean -xu > .c
diff .a .b
diff .b .c


# same for h3
flutile trim h3-ha1 --conversion dna2aa  3.2010.2.fna 2> /dev/null | smof clean -xu > .a
flutile trim h3-ha1 --conversion dna2dna 3.2010.2.fna 2> /dev/null | smof translate | smof clean -xu > .b
flutile trim h3-ha1 --conversion aa2aa   3.2010.2.faa 2> /dev/null | smof clean -xu > .c
diff .a .b
diff .b .c

rm -f .a .b .c
