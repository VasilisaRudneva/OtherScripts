#!/bin/bash

if [ $# -ne 2 ]
then
    echo "Usage: <amp> <del>"
    exit -1
fi
a=$1
d=$2

echo "amp_formula='${a}'; del_formula='${d}'" | cat - PrepareROCInput.R | R --vanilla --slave &
