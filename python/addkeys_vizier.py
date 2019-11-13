# add keywords to the file to make it compatible with vizier


import os
import glob
import astropy
from astropy.io import fits

def addkeys_agns(filename):

    hdu = fits.open(filename)
    header = hdu[1].header
    Ncols = int(hdu[1].header['TFIELDS']) 
    Nfreq=Ncols-12
    
    header['TDISP1'] = 'F6.4'
    #header['TUNIT1'] = 'log (erg/s/Hz)'
    for i in range(1, 1+Nfreq):
        header['TDISP%i'%(i+1)] = 'F30.15'
        #header['TUNIT%i'%(i+1)] = 'mJy'
    i=i+1
    header['TDISP%i'%(i+1)] = 'F6.4' #Mh
    #header['TUNIT%i'%(i+1)] = 'log(Msun)' #Mh
    i=i+1
    header['TDISP%i'%(i+1)] = 'F8.4' #x_coord
    #header['TUNIT%i'%(i+1)] = 'degs' 
    i=i+1
    header['TDISP%i'%(i+1)] = 'F8.4' #y_coord
    #header['TUNIT%i'%(i+1)] = 'degs' 
    i=i+1
    header['TDISP%i'%(i+1)] = 'F8.4' #latitude

    #header['TUNIT%i'%(i+1)] = 'degs' 
    i=i+1
    header['TDISP%i'%(i+1)] = 'F8.4' #longitude
    #header['TUNIT%i'%(i+1)] = 'degs' 
    i=i+1
    header['TDISP%i'%(i+1)] = 'F8.7' #redshift
    #header['TUNIT%i'%(i+1)] = 'none' 
    i=i+1
    header['TDISP%i'%(i+1)] = 'F15.10' #phys size
    #header['TUNIT%i'%(i+1)] = 'Kpc' 
    i=i+1
    header['TDISP%i'%(i+1)] = 'F15.10' #angle
    #header['TUNIT%i'%(i+1)] = 'degs' 
    i=i+1
    header['TDISP%i'%(i+1)] = 'F15.10' #size
    #header['TUNIT%i'%(i+1)] = 'arcsec' 
    i=i+1
    header['TDISP%i'%(i+1)] = 'F8.7' #Rs
    #header['TUNIT%i'%(i+1)] = 'none' 
    i=i+1
    header['TDISP%i'%(i+1)] = 'F1.0' #popflag
    #header['TUNIT%i'%(i+1)] = 'none' 
    hdu.writeto(filename,overwrite=True)
    return



def addkeys_sfgs(filename):
   

    print(filename)
    hdu = fits.open(filename)
    header = hdu[1].header
    Ncols = int(hdu[1].header['TFIELDS']) 
    Nfreq=Ncols-11
    print(Nfreq)
    header['TDISP1'] = 'F6.4'
    #header['TUNIT1'] = 'log (Msun/yr)'
    for i in range(1, 1+Nfreq):
        header['TDISP%i'%(i+1)] = 'F30.15'
    #    header['TUNIT%i'%(i+1)] = 'mJy'
    i=i+1
    header['TDISP%i'%(i+1)] = 'F6.4' #Mh
    #header['TUNIT%i'%(i+1)] = 'log(Msun)' #Mh
    i=i+1
    header['TDISP%i'%(i+1)] = 'F8.4' #x_coord
    #header['TUNIT%i'%(i+1)] = 'degs' 
    i=i+1
    header['TDISP%i'%(i+1)] = 'F8.4' #y_coord
    #header['TUNIT%i'%(i+1)] = 'degs' 
    i=i+1
    header['TDISP%i'%(i+1)] = 'F8.4' #latitude
    #header['TUNIT%i'%(i+1)] = 'degs' 
    i=i+1
    header['TDISP%i'%(i+1)] = 'F8.4' #longitude
    #header['TUNIT%i'%(i+1)] = 'degs' 
    i=i+1
    header['TDISP%i'%(i+1)] = 'F8.7' #redshift
    #header['TUNIT%i'%(i+1)] = 'none' 
    i=i+1
    header['TDISP%i'%(i+1)] = 'F15.10' #size
    #header['TUNIT%i'%(i+1)] = 'arcsec' #size
    i=i+1
    header['TDISP%i'%(i+1)] = 'F8.7' #e1
    #header['TUNIT%i'%(i+1)] = 'none' 
    i=i+1
    header['TDISP%i'%(i+1)] = 'F8.7' #e2
    #header['TUNIT%i'%(i+1)] = 'none' 
    i=i+1
    header['TDISP%i'%(i+1)] = 'F1.0' #popflag
    #header['TUNIT%i'%(i+1)] = 'none' 

    hdu.writeto(filename,overwrite=True)


    return

filename='/home/a.bonaldi/local2/scratch/Bonaldi/Radio_srccnt/runs_delivery/Sept2019/catalogue_SFGs_complete_medium.fits'
addkeys_sfgs(filename)

filename='/home/a.bonaldi/local2/scratch/Bonaldi/Radio_srccnt/runs_delivery/Sept2019/catalogue_AGNs_complete_deep.fits'
#addkeys_agns(filename)
