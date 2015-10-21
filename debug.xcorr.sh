#!/bin/bash

## compute spatial correlation function
##  with TPACF code (Dolence & Brunner 2007)
## b. wibking, jun 2015

TPACF=$HOME/tpacf

DOSPATIAL=1
JACKKNIFE_SAMPLES=16

DOAUTOCORR=0
DOCROSSCORR=1

NBINS=30
MINBIN=0.1
MAXBIN=30
BOXSIZE=1911.2

THREADS=16
TEMPFILE=sfilelist.xcorr
TEMPFILE2=sparams.xcorr

echo "$1 1
$2 1" > $TEMPFILE

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

echo "Converting HDF5 to unformatted binary"
gdb --args $TPACF/bin/precompute $TEMPFILE $DOSPATIAL $JACKKNIFE_SAMPLES

echo "Computing spatial correlation"
OMP_NUM_THREADS=$THREADS gdb --args $TPACF/bin/correlate $TEMPFILE2



