#!/usr/bin/env bash

set -e
set -u

diff <(flutile aadiff --subtype=H1 h1-aadiff-ex.faa) h1-aadiff-ex.txt

flutile aadiff --caton82 --annotation-tables=rbs.txt --subtype=H1 --join-annotations 1B.2.1-aadiff-test.faa > .obs-aadiff
diff .obs-aadiff .exp-aadiff
