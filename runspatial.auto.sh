#!/bin/bash

## compute spatial correlation function
##  with TPACF code (Dolence & Brunner 2007)
## b. wibking, jun 2015

TPACF=$HOME/tpacf

DOSPATIAL=1
JACKKNIFE_SAMPLES=16

DOAUTOCORR=1
DOCROSSCORR=0

NBINS=30
MINBIN=0.1
MAXBIN=30
BOXSIZE=1911.2

THREADS=16
TEMPFILE=sfilelist.xcorr
TEMPFILE2=sparams.xcorr

echo "$1 1
$1 1" > $TEMPFILE

echo "$TEMPFILE
$DOAUTOCORR
$DOCROSSCORR
0
5
$NBINS
$MINBIN
$MAXBIN
1
$BOXSIZE" > $TEMPFILE2

echo "Computing correlation..."
OMP_NUM_THREADS=$THREADS $TPACF/bin/correlate $TEMPFILE $TEMPFILE2 $DOSPATIAL $JACKKNIFE_SAMPLES



