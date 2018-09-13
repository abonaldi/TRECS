# list all output files produced by sampler on a given folder
# prepare input file for wrapper
# collate all files into one catalogue for AGNs and one for SFGs
# add keyword to header to make output compatible with Vizier

import os
import glob
import astropy

def collate_sfgs(folder1,folder2,lat_target,lon_target,sim_side,do_clustering,tag):
    

    from astropy.io import fits
    file=folder2+'/catalogue_SFGs_complete_'+tag+'.fits'

    try:
        os.system('rm '+file)
    except OSError:
        pass

    results0=glob.glob(folder1+'/*UVgal.fits')
    results1=glob.glob(folder1+'/*spheroids.fits')
    results2=glob.glob(folder1+'/*lens_spher.fits')

    results=results0+results1+results2


    hdul=fits.open(results[0])
    Ncols = int(hdul[1].header['TFIELDS']) 
    nfiles=len(results)

    print('Ncols',Ncols)
    print('Nfiles',nfiles)

    input_file = open('infile_wrapper_SFGs_'+tag+'.ini','w')
    input_file.write('Nfiles=%i\n'%(nfiles))
    input_file.write('Ncols=%i\n'%(Ncols))
    input_file.write('lat_target=%f\n'%(lat_target))
    input_file.write('lon_target=%f\n'%(lon_target))
    input_file.write('sim_side=%f\n'%(sim_side))
    input_file.write('do_clustering=%i\n'%(do_clustering))

    for i in range(nfiles):
        input_file.write('cat%i=%s\n'%(i+1,results[i]))

    input_file.write('outcat=%s\n'%file)
    input_file.close()

    return



def collate_agns(folder1,folder2,lat_target,lon_target,sim_side,do_clustering,tag):

    from astropy.io import fits
    file=folder2+'/catalogue_AGNs_complete_'+tag+'.fits'
    
    try:
        os.system('rm '+file)
    except OSError:
        pass

    results0=glob.glob(folder1+'/*SS_AGN.fits')
    results1=glob.glob(folder1+'/*FSRQ.fits')
    results2=glob.glob(folder1+'/*BLLac.fits')
    
    results=results0+results1+results2

    print(results[0])

    hdul=fits.open(results[0])
    Ncols = int(hdul[1].header['TFIELDS'])   
    
    nfiles=len(results)

    input_file = open('infile_wrapper_AGNs_'+tag+'.ini','w')
    input_file.write('Nfiles=%i\n'%(nfiles))
    input_file.write('Ncols=%i\n'%(Ncols))
    input_file.write('lat_target=%f\n'%(lat_target))
    input_file.write('lon_target=%f\n'%(lon_target))
    input_file.write('sim_side=%f\n'%(sim_side))
    input_file.write('do_clustering=%i\n'%(do_clustering))

    for i in range(nfiles):
        input_file.write('cat%i=%s\n'%(i+1,results[i]))
    input_file.write('outcat=%s\n'%file)
    input_file.close()


    return

folder2='/home/a.bonaldi/data_challenges/inputs/T-RECS_cats/'
folder1='/home/a.bonaldi/data_challenges/inputs/T-RECS_cats/run1/'
tag='chall_1'

lat_target=-30.
lon_target=0.

do_clustering=1
sim_side=3.

#folder2='../'
#folder1='../tests/'
#tag='test'

#lat_target=0.
#lon_target=0.

#do_clustering=1
#sim_side=5.


try:
    os.mkdir(folder2)
except OSError:
    pass



collate_agns(folder1,folder2,lat_target,lon_target,sim_side,do_clustering,tag)
collate_sfgs(folder1,folder2,lat_target,lon_target,sim_side,do_clustering,tag)

