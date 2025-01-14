#!/bin/bash
if [ -z "$1" ]; then
    echo "Usage: $0 <variable_name>"
    exit 1
fi
gnuplot -e "variable_name='$1'" plot.gp
