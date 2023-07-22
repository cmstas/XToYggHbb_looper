#!/bin/bash

DIR=$1
YEAR=$2 # all, 2018, 2017, 2016APV, 2016nonAPV
DATA=$3 # 0, 1
BKG=$4 # 0, 1
SIG=$5 # 0, 1
SAM=$6
LMM=$7 # 0, 1
while ! [ -z "$8" ]; do
    FLAGS="$FLAGS $8"; shift;
done
COMMAND="./main.exe $DIR $YEAR $DATA $BKG $SIG $SAM $LMM$FLAGS"

echo ""
echo "Arguments:"
echo "----------"
echo "Directory = "$DIR
echo "Year = "$YEAR
echo "Run data = "$DATA
echo "Run bkg = "$BKG
echo "Run signal = "$SIG
echo "Sample = "$SAM
echo "Low mass mode = "$LMM
echo "Variation flags = "${FLAGS#", "}
echo ""

eval $COMMAND
