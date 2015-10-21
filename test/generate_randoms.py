#!/usr/bin/env python

## Generate a spatial Poisson distribution of points for
## autocorrelation and crosscorrelation tests

import numpy as np
import pandas as pd
from odo import odo

def randomPointsInBox(boxLength,npoints):
    randomArray = np.random.uniform(low=0.0,high=boxLength,size=(npoints,3))
    return pd.DataFrame(randomArray,columns=['x','y','z'])

def singlePointInBox(boxLength):
    randomArray = np.array([[boxLength/2.,boxLength/2.,boxLength/2.],])
    print randomArray.shape
    return pd.DataFrame(randomArray,columns=['x','y','z'])

def main():
    import argparse
    parser = argparse.ArgumentParser(description='generate hdf5 file of random points.')
    parser.add_argument('boxLength',type=float)
    parser.add_argument('npoints',type=int)
    parser.add_argument('outputfile')
    args = parser.parse_args()

    if args.npoints>1:
        randomDf = randomPointsInBox(args.boxLength,args.npoints)
    else:
        randomDf = singlePointInBox(args.boxLength)

    fileOutput = odo(randomDf,args.outputfile+"::/particles")
    print fileOutput

main()
