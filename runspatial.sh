#!/bin/bash
# compute spatial correlation function
DOSPATIAL=1
JACKKNIFE_SAMPLES=0
THREADS=16

echo
bin/precompute sfilelist $DOSPATIAL $JACKKNIFE_SAMPLES
echo "Computing spatial correlation"
#TODO better describe parameters in sparams.in
OMP_NUM_THREADS=$THREADS bin/correlate sparams.in



