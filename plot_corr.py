#import pyqtgraph as plt
import matplotlib.pyplot as plt
import numpy as np
import math

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('hod_filename', help='filename of output DD pair counts file')
parser.add_argument('dm_filename', help='filename of output DD pair counts file')
parser.add_argument('xcorr_filename', help='filename of output cross-correlation pair counts file')

args = parser.parse_args()

filenames = [(args.hod_filename,1), (args.dm_filename,1), (args.xcorr_filename,2)]
x = []
y = []

for filename,skiprows in filenames:
	f = open(filename,'r')
	header = f.readline()
	#f.close()

	header_parts = header.split(' ')
	NumBins = header_parts[0]
	NumJackknife = header_parts[1]
	npart = header_parts[2]
	NumBins = int(NumBins) # ~30
	npart = int(npart) # 10^6

	print NumBins,NumJackknife,npart
	
	if skiprows>1:
		header_2 = f.readline()
		header2_parts = header_2.split(' ')
		NumPart2 = header2_parts[0]
		print NumPart2
	f.close()

	part = np.loadtxt(filename,skiprows=skiprows)

	maxDist = 30.0 # can this be determined from the file?
	minDist = 0.1

	Lbox = 1911.2 # TODO: add as an argument

	bins = [maxDist*maxDist*(minDist/maxDist)**((2.*i)/(NumBins)) for i in xrange(NumBins)]
	bins.append(minDist*minDist)
	rbins = np.array(bins)**0.5
	rbins = rbins[::-1]
	part = part[::-1]

	if skiprows==1:
		ndens = npart / (Lbox**3)
	else:
		ndens = int(NumPart2)/(Lbox**3)
	
	V = rbins[1:]**3 - rbins[:-1]**3
	RR = 4./3. * math.pi * V *ndens*npart
	allxsi = part[:,1]/RR - 1.0

	x.append(part[:,0])
	y.append(allxsi)

	#plt.plot(np.log10(part[:,0]),np.log10(allxsi),'-o')
	plt.plot(part[:,0],allxsi,'-o',label=filename)
	plt.xscale('log')
	plt.yscale('log')

bins = x[0]
xsi = y[0]
xsi_DM = y[1]
xsi_xcorr = y[2]	

plt.xlim((bins[0],bins[-1]))
plt.legend(loc='best')

# compute galaxy bias                                                                                
bias = np.sqrt(xsi/xsi_DM)

# compute galaxy-mass 'correlation coefficient'                                                      
corr = xsi_xcorr/np.sqrt(xsi*xsi_DM) # can be (much) greater than 1 (but this is still correct)!	

plt.figure()
plt.plot(bins, bias, '-o', label="bias")
plt.plot(bins, corr, '-o', label="pseudo-correlation")
plt.xscale('log')
plt.xlim((bins[0],bins[-1]))

plt.legend(loc='best')
plt.show()
