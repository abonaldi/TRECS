from sklearn.neighbors import NearestNeighbors
import numpy as np
import os, glob, sys
from astropy.io import fits
import matplotlib
from matplotlib import pyplot as plt
import astropy
from collections import Counter
from astropy.cosmology import LambdaCDM
import time
 
tstart = time.time()
 
 
 
cat_name1 = '/home/a.bonaldi/local2/scratch/Bonaldi/Radio_srccnt/TESTS_newtrecs/SDC2/test2/hi/catalogue_HI_z0.10.fits'
cat_fits1 = fits.open(cat_name1)
cat1 = cat_fits1[1].data
cols1 = cat_fits1[1].columns.names


cat_name2 = '/home/a.bonaldi/local2/scratch/Bonaldi/Radio_srccnt/TESTS_newtrecs/SDC2/test2/continuum/catalogue_continuum_z0.10.fits'
cat_fits2 = fits.open(cat_name2)
cat2 = cat_fits2[1].data
cols2 = cat_fits2[1].columns.names
 
print ((cols1))
 


MHI = cat1['MHI']

MHI_2 = cat2['MHI_pred']

attr1=np.array([MHI]).T
attr2=np.array([MHI_2]).T

print ('this is the shape')
print (attr1)

nbrs = NearestNeighbors(n_neighbors=1, algorithm='kd_tree').fit(attr1)
distances, indices = nbrs.kneighbors(attr2)

print (indices.shape)
print (attr1.shape)
print (attr2.shape)

print ('done')
 
 
#attr1 = np.vstack((MHI, MH, r_HI, incl)).T  #collate different vecorts in an array and transpose because we want them in column order

#attr2 = np.vstack((MHI, MH, r_HI, incl)).T

 
 
 
 
#print (attr1.shape)
 
 
 
 
 
#nbrs = NearestNeighbors(n_neighbors=1, algorithm='kd_tree').fit(attr1)
#distances, indices = nbrs.kneighbors(attr2)


print (distances)
print (indices)
 
tend = time.time()
print ('...done in {0} seconds.'.format(tend-tstart))
