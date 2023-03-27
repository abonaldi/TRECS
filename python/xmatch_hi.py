
# Assigns counterparts between HI and continuum T-RECS catalogues by creating a crossmatched catalogue
# the quantity matched between the two catalogues if the HI mass, called MHI in the HI catalogue and MHI_pred in the continuum catalogue
# the output is a catalogue containing the fields of both catalogues, with matched objects on the same row
# if one of the two catalogues in not present, a crossmatched catalogue is produced with the same data as the original ones but with the same format as the crossmatched ones (extra columns), with a -100 flag for the fields that are not present

# History
#---------------------------------
# A. Bonaldi 14/7/21 first version





from sklearn.neighbors import NearestNeighbors
import numpy as np
import os, glob, sys
from astropy.io import fits
#import matplotlib
#from matplotlib import pyplot as plt
import astropy
#from collections import Counter
from astropy.table import Table
from astropy.table import Column

import time
 
tstart = time.time()
print('Number of arguments:', len(sys.argv), 'arguments.')
print('Argument List:', sys.argv)
zmin=float(sys.argv[1])
zmax=float(sys.argv[2])
print(zmin,zmax)

fov=5. # Field of view of the crossmatched catalogue

path_HI='/home/a.bonaldi/data-cold-for-backup/Radio_srccnt/runs_paper2/hi_rev/' #path if HI catalogues
tag='HI'      #tag in the file name of HI catalogues
sim_side=5.   #size of square HI simulation (degs)

path_cont='/home/a.bonaldi/data-cold-for-backup/Radio_srccnt/runs_paper2/continuum/' #path of continuum catalogue
sim_side_c=5. #size of square continuum simulation (degs)
tag_c='continuum' #tag in the file name of continuum catalogues

path_out='/home/a.bonaldi/data-cold-for-backup/Radio_srccnt/runs_paper2/cross_rev/' #path of output catalogues

####end general settings

if (sim_side < fov) or (sim_side_c < fov):
    print("Error: input catalogue(s) smaller than the required FoV")
    exit()

halfside=fov/2.


redshift_names=['0.01','0.02','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00','1.20','1.40','1.60','1.80','2.00','2.20','2.40','2.60','2.80','3.00','3.20','3.40','3.60','3.80','4.00','4.20','4.40','4.60','4.80','5.00','5.20','5.40','5.60','5.80','6.00','6.20','6.40','6.60','6.80','7.00','7.20','7.40','7.60','7.80','8.00','8.20','8.40','8.60','8.80','9.00','9.20','9.40','9.60','9.80','10.0']


# initialise format
#read in file headers to format output correctly

for i in range(len(redshift_names)):
    z=redshift_names[i]
    cat_name1 = path_HI+'catalogue_'+tag+'_z'+z+'.fits'
    if (os.path.isfile(cat_name1) == True):
        cat_fits1 = fits.open(cat_name1)
        cols1 = cat_fits1[1].columns.names
        break

for i in range(len(redshift_names)):
    z=redshift_names[i]
    cat_name2 = path_cont+'catalogue_'+tag_c+'_z'+z+'.fits'
    if (os.path.isfile(cat_name2) == True):
        cat_fits2 = fits.open(cat_name2)
        cols2 = cat_fits2[1].columns.names
        break

#modify the name of columns that are present in both catalogues



for i in range(len(cols2)):
    if cols2[i] in cols1:
        cols2[i] = cols2[i]+'_1'


#end initialise format

#start main loop
for i in range(len(redshift_names)):


    z=redshift_names[i]

    if (np.float(z) >= zmin) and (np.float(z) <= zmax):

        print('********************')
        print('Processing redshift',z)
        print('********************')

        ngals=0 #initialise number of objects - for the case where file does not exist
        nhaloes=0
        
        cat_name1 = path_HI+'catalogue_'+tag+'_z'+z+'.fits'
        
        if (os.path.isfile(cat_name1) == True):
            cat1 = Table.read(cat_name1, format='fits')
        
            #select FoV only
            cat1=cat1[(np.abs(cat1['x_coord']) <= halfside)*(np.abs(cat1['y_coord']) <= halfside)]     
            ngals=len(cat1['x_coord'])
        
        cat_name2 = path_cont+'catalogue_'+tag_c+'_z'+z+'.fits'
        
        if (os.path.isfile(cat_name2) == True):
            cat2 = Table.read(cat_name2, format='fits')
            cat_fits2 = fits.open(cat_name2)
            cols2_old = cat_fits2[1].columns.names

                        
            #rename columns
            for col in range(len(cols2)):
                cat2.rename_column(cols2_old[col] , cols2[col])
            
            #here select only the portion of the catalogue that we need
            cat2=cat2[(np.abs(cat2['x_coord_1']) <= halfside)*(np.abs(cat2['y_coord_1']) <= halfside)] 
            nhaloes=len(cat2['x_coord_1'])


        
        # bring catalogues to a common format by adding missing columns with a -100 flag value
        if (ngals !=0):
            emptycol=np.zeros(ngals)-100
            for i, col in enumerate(cols2):
                cat1.add_column(Column(name=cols2[i],data=emptycol))


        
        if (nhaloes !=0):
            emptycol=np.zeros(nhaloes)-100
        
            for i, col in enumerate(cols1):
                cat2.add_column(Column(name=cols1[i],data=emptycol),index=i)
            # fill lat, lon, redshift with fields _1 if primary is missing
            cat2['x_coord']=cat2['x_coord_1']
            cat2['y_coord']=cat2['y_coord_1']
            cat2['redshift']=cat2['redshift_1']
            cat2['Mh']=cat2['Mh_1']
        

        cat_name1_out = path_out+'catalogue_HI_continuum_z'+z+'.fits' 
       

        #start reading catalogue values

        if (ngals !=0) and (nhaloes !=0.):
            M1 = cat1['MHI']
            M2 = cat2['MHI_pred']

            minmass=np.min(M1)-1. #this is the minimum mass to be matched, decreased to be conservative
            
            attr1=np.array([M1]).T   #observed catalogue
            attr2=np.array([M2]).T   #lightcone

        
            print('number of HI galaxies',ngals)
            print('number of continuum galaxies',nhaloes)

       

            ###### ANALYSIS STARTS
            print('Start nearest neighbour')
            nsample=np.min([20,nhaloes]) #number of nearest neighbour considered for matching
            nbrs = NearestNeighbors(n_neighbors=nsample, algorithm='kd_tree').fit(attr2)
            distances, indices = nbrs.kneighbors(attr1)
            print('Neaerst neighbour done')
        
            print('number of matches',indices.shape)
        

            #sort to match high mass first - more rare and difficult to match
            # go trough suggested matches and assign 
        
            indices_sortM1=np.flip(np.argsort(M1)) # indices for observed masses from largest to smallest
        

            for k in range(ngals):

                ###this is the galaxy we want to associate with an halo
                j=indices_sortM1[k]

                for jj in range(nsample):
                        
                    if (M2[indices[j,jj]] != -100.):
                        
                        for i, col in enumerate(cols2):
                            
                            cat1[cols2[i]][j] = cat2[cols2[i]][indices[j,jj]]  #copy this source to the spare columns of the counterpart in catalogue 1
                            cat2[cols2[i]][indices[j,jj]]=-100. #flag this source as matched in catalogue 2
                                                                                 
                        jj_match=jj  #this is to flag the mass later to avoid repetitions
                      
                        
                        break
                
                        M2[indices[j,jj_match]]=-100. #flagging this mass so it cannot be reused
                


            #append non-matched cat2 source at the end of cat1                 
            empty=cat2['x_coord_1'] #column used to check if this galaxy has been matched 
            for i in range(nhaloes):
                if (empty[i] != -100.):
                    row=cat2[i]
                    cat1.add_row(row)
                    
        if (ngals+nhaloes !=0):
            catout=Table()
            
            if (ngals ==0) and (nhaloes != 0):
                catout=cat2 #rewrite cat2, but with the extra columns 
            else:
                catout=cat1    

            # remove duplicate coordinates and redshifts from the final catalogue.
            # Those were all random uniform generation within the given ranges so
            # the redundant set do not contain any useful information
            # the second dark halo mass Mh_1 is kept as it comes from continuum model instead of HI 
            catout.remove_columns(['x_coord_1','y_coord_1','redshift_1','latitude_1','longitude_1'])
                      
            print('writing updated catalogue file')
            catout.write(cat_name1_out,format='fits', overwrite = True)

                
tend = time.time()
print ('...done in {0} seconds.'.format(tend-tstart))
