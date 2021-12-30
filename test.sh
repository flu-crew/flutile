#!/usr/bin/env bash

set -e
set -u
set -v

cd test-data
./test.sh

cd clades
./test.sh

cd ../ha1
./test.sh
