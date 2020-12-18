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
