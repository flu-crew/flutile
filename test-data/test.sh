#!/usr/bin/env bash

set -e
set -u

diff <(flutile aadiff --subtype=H1 h1-aadiff-ex.faa 2> /dev/null) h1-aadiff-ex.txt
