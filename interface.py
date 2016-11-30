#!/usr/bin/python
# -*- coding: utf-8 -*-
from astropy.io import fits
from sgpdecon import *
import numpy as np 
import glob,os


########################################################################
"""This is the interface part of fits files deconvolution program, set 
correct path as "filepath/psf" and proper value as iteration times. """
########################################################################
"Set Default"
path="G:"
path2="G:"
filepath=path+"/*.fit*"
decpath=path+"/deconvolution/"
os.popen('rm -rf '+decpath)
os.mkdir(decpath)
iteration=raw_input("Iteration= ")
iteration=int(iteration)

if __name__ =='__main__':
	fitlist=glob.glob(filepath)
	declist=fitlist[:]
for i in range(1,len(fitlist)+1):
	print i,fitlist[i-1]
	hdu=fits.open(fitlist[i-1])
	gray=hdu[0].data
	h=hdu[0].header['NAXIS1']
	w=hdu[0].header['NAXIS2']
	psfload=path2+"/"+"fivehundred11times11psf"
	psfdata=np.loadtxt(psfload)
	datetype=type(psfdata[0,0])
	psf=np.zeros([h,w],dtype=datetype)
	psf[505:516,505:516]=psfdata
	#gray=np.transpose(gray)
	if __name__ =='__main__':
		parameters=[gray,psf,iteration]
		newimg=sgpdecon(parameters)
		dec = fits.PrimaryHDU(newimg)
		decpic = fits.HDUList([dec])
		decsave=decpath+declist[i-1].split("/")[-1]
		decpic.writeto(decsave)
		

