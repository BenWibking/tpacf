import matplotlib.pyplot as plt
import numpy as np
import math

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('redshift',type=float)
parser.add_argument('logMmin',type=float)
parser.add_argument('siglogM',type=float)
parser.add_argument('logM1',type=float)
parser.add_argument('logM0',type=float)
parser.add_argument('alpha',type=float)
parser.add_argument('mindist',type=float)
parser.add_argument('maxdist',type=float)
parser.add_argument('boxsize',type=float)
parser.add_argument('hod_filename', help='filename of output DD pair counts file')
parser.add_argument('dm_filename', help='filename of output DD pair counts file')
parser.add_argument('xcorr_filename', help='filename of output cross-correlation pair counts file')

args = parser.parse_args()

filenames = [(args.hod_filename,1,'galaxy'), (args.dm_filename,1,'matter'), (args.xcorr_filename,2,'galaxy-matter')]
x = []
y = []

fig = plt.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8])

for filename,skiprows,title in filenames:
	f = open(filename,'r')
	header = f.readline()

	header_parts = header.split(' ')
	NumBins = header_parts[0]
	NumJackknife = header_parts[1]
	npart = header_parts[2]
	NumBins = int(NumBins)
	npart = int(npart)

	print "numbins:",NumBins,"numjackknife:",NumJackknife,"npart:",npart
	
	if skiprows>1:
		header_2 = f.readline()
		header2_parts = header_2.split(' ')
		NumPart2 = header2_parts[0]
		print "Numpart2:",NumPart2
	f.close()

	part = np.loadtxt(filename,skiprows=skiprows)

	#maxDist = 30.0 # can this be determined from the file?
	#minDist = 0.1

	minDist = args.mindist
	maxDist = args.maxdist
	Lbox = args.boxsize

	bins = [maxDist*(minDist/maxDist)**((1.*i)/(NumBins)) for i in xrange(NumBins)]
	bins.append(minDist)
	rbins = np.array(bins)
	rbins = rbins[::-1]
	part = part[::-1]

	if skiprows==1:
		ndens = npart / (Lbox**3)
		mylabel = r"{t} $n={ndens:.2e}$".format(t=title,ndens=ndens)
	else:
		ndens = int(NumPart2)/(Lbox**3)
		mylabel = r"{t}".format(t=title)

	V = rbins[1:]**3 - rbins[:-1]**3
	RR = 4./3. * math.pi * V *ndens*npart
	allxsi = part[:,1]/RR - 1.0

	x.append(part[:,0])
	y.append(allxsi)

	ax.plot(part[:,0],allxsi,'-o',label=mylabel)
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlabel(r'r ($h^{-1}$ Mpc)')

bins = x[0]
xsi = y[0]
xsi_DM = y[1]
xsi_xcorr = y[2]	

## save the properly-normalized values to HDF5 file
from odo import odo
output_filename = args.hod_filename+".2pcf.hdf5"

odo(bins,output_filename+'::/bins').file.close()
odo(xsi,output_filename+'::/xsi_gal').file.close()
odo(xsi_DM,output_filename+'::/xsi_matter').file.close()
odo(xsi_xcorr,output_filename+'::/xsi_gal_matter').file.close()

ax.set_xlim((bins[0],bins[-1]))
ax.legend(loc='best')
hod_label = r"$\log M_{{min}}={logMmin}$, $\sigma_{{\log M}}={siglogM}$, $\log M_1={logM1}$, $\log M_0={logM0}$, $\alpha={alpha}$".format(siglogM=0.5, logMmin=12.5, logM0=12.5, logM1=13.5, alpha=1.0)
ax.text(.02,.01,hod_label,verticalalignment='bottom',horizontalalignment='left',fontsize=16,transform=ax.transAxes)
sim_label = r"$z={z}$ Abacus simulation".format(z=0.57)
ax.text(.02,.99,sim_label,verticalalignment='bottom',horizontalalignment='left',fontsize=16,transform=ax.transAxes)
plt.suptitle(r'two-point real-space correlations')

fig.savefig('2pcf_gal.pdf')

# compute galaxy bias                                                                                
bias = np.sqrt(xsi/xsi_DM)

# compute galaxy-mass 'correlation coefficient'                                                      
corr = xsi_xcorr/np.sqrt(xsi*xsi_DM) # can be (much) greater than 1 (but this is still correct)!	

## save to HDF5 file
odo(bias,output_filename+'::/galaxy_bias').file.close()
odo(corr,output_filename+'::/gal_matter_correlation').file.close()

import h5py
with h5py.File(output_filename) as h5f:
	h5f['bins'].attrs.create('redshift',args.redshift)
	h5f['bins'].attrs.create('logMmin',args.logMmin)
	h5f['bins'].attrs.create('logM0',args.logM0)
	h5f['bins'].attrs.create('logM1',args.logM1)
	h5f['bins'].attrs.create('siglogM',args.siglogM)
	h5f['bins'].attrs.create('alpha',args.alpha)

## plotting

plt.figure()
plt.plot(bins, bias, '-o', label="galaxy bias")
plt.xscale('log')
plt.xlim((bins[0],bins[-1]))
plt.xlabel(r'r ($h^{-1}$ Mpc)')
plt.legend(loc='best')
plt.savefig('bias_gal.pdf')

plt.figure()
plt.plot(bins, corr, '-o', label="galaxy-matter pseudo-correlation")
plt.xscale('log')
plt.xlabel(r'r ($h^{-1}$ Mpc)')
plt.xlim((bins[0],bins[-1]))
plt.legend(loc='best')
plt.savefig('r_gm.pdf')

plt.show()
