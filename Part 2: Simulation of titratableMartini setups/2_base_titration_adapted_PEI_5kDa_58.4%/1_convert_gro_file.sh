#!/bin/bash

run="python3 ../scripts/convert_gro_to_titratable.py"

$run -f solvated.gro -o temp.gro -sel "name TN6d or name SN4 or name SN3a" -bead base 
$run -f temp.gro -o start.gro -sel "name W" -bead water

rm -f temp.gro
