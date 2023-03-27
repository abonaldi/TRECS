# list all output files produced by sampler on a given folder
# prepare input file for wrapper
# collate all files into one catalogue for AGNs and one for SFGs
# add keyword to header to make output compatible with Vizier

import os
import glob
import astropy


def collate(folder1,folder2,lat_target,lon_target,sim_side,do_clustering,tag):
    

    from astropy.io import fits
    file=folder2+'/catalogue_complete_'+tag+'.fits'

    try:
        os.system('rm '+file)
    except OSError:
        pass

    results0=glob.glob(folder1+'/*.fits')

    results=results0
    

    hdul=fits.open(results[0])
    Ncols = int(hdul[1].header['TFIELDS']) 
    nfiles=len(results)

    print('Ncols',Ncols)
    print('Nfiles',nfiles)

    input_file = open('infile_wrapper_'+tag+'.ini','w')
    input_file.write('Nfiles=%i\n'%(nfiles))
    input_file.write('Ncols=%i\n'%(Ncols))
    input_file.write('lat_target=%f\n'%(lat_target))
    input_file.write('lon_target=%f\n'%(lon_target))
    input_file.write('sim_side=%f\n'%(sim_side))
    input_file.write('do_clustering=%i\n'%(do_clustering))

    for i in range(nfiles):
        input_file.write('cat%i=%s\n'%(i+1,results[i]))

    input_file.write('outcat=%s\n'%file)
    input_xfile.close()

    return




#folder1='/home/a.bonaldi/data-cold-for-backup/Radio_srccnt/runs_paper2/continuum_wide/'
#folder2='/home/a.bonaldi/data-cold-for-backup/Radio_srccnt/runs_paper2/'

#tag='continuum_wide_10Mar'

#folder1='/home/a.bonaldi/data-cold-for-backup/Radio_srccnt/runs_paper2/continuum_theta_p/'
#folder2='/home/a.bonaldi/data-cold-for-backup/Radio_srccnt/runs_paper2/'

#tag='continuum_thetap'



folder1='/home/a.bonaldi/data-cold-for-backup/Radio_srccnt/runs_paper2/cross_rev/'
folder2='/home/a.bonaldi/data-cold-for-backup/Radio_srccnt/runs_paper2/'

tag='cross_rev'

#folder1='/home/a.bonaldi/data-cold-for-backup/Radio_srccnt/runs_paper2/cross_test/'
#folder2='/home/a.bonaldi/data-cold-for-backup/Radio_srccnt/runs_paper2/'

#tag='cross_6dec'


#folder1='/home/a.bonaldi/data-cold-for-backup/SDC3/TRECS_cats/run3/8/'
#folder2='/home/a.bonaldi/data-cold-for-backup/SDC3/TRECS_cats/'

#tag='SDC3_v3_8'

#folder1='/home/a.bonaldi/data-cold-for-backup/Radio_srccnt/runs_magnetism/run1/'
#folder2='/home/a.bonaldi/data-cold-for-backup/Radio_srccnt/runs_magnetism/'



#tag='magnetism_v1'


#folder1='/home/a.bonaldi/data-cold-for-backup/Radio_srccnt/runs_lsync/test1/'
#folder2='/home/a.bonaldi/data-cold-for-backup/Radio_srccnt/runs_lsync/'

#tag='test1'


#tag='100nJy'
#folder1='/home/a.bonaldi/data-cold-for-backup/Radio_srccnt/runs_sensitivity_calculator/limit_100/'
#folder2='/home/a.bonaldi/data-cold-for-backup/Radio_srccnt/runs_sensitivity_calculator/'


lat_target=-30.
lon_target=0.

do_clustering=1
#sim_side=25
sim_side=5


try:
    os.mkdir(folder2)
except OSError:
    pass



collate(folder1,folder2,lat_target,lon_target,sim_side,do_clustering,tag)

