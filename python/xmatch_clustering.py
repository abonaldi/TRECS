# 16/9/20
# overwrites the sky coordinates (x_coord, y_coord, redshift) of galaxies with those of a lightcone of a cosmological simulation, thus reproducing clustering
# matching galaxies between the DM catalogue and the observed catalogue based on DM halo mass
# use nearest neighbours for associating mass
# option to use local density as well to introduce environmental dependencies. 




from sklearn.neighbors import NearestNeighbors
import numpy as np
import os, glob, sys
from astropy.io import fits
import matplotlib
from matplotlib import pyplot as plt
import astropy
from collections import Counter
#from astropy.cosmology import LambdaCDM
from astropy.table import Table
#from astropy.table import hstack
#from astropy.table import vstack

import time
 
tstart = time.time()
print('Number of arguments:', len(sys.argv), 'arguments.')
print('Argument List:', sys.argv)
#exit() 
###general settings
#zmin=7.
#zmax=7.2
zmin=float(sys.argv[1])
zmax=float(sys.argv[2])
print(zmin,zmax)

#path='/home/a.bonaldi/local2/scratch/Bonaldi/Radio_srccnt/TESTS_newtrecs/SDC2/test4/continuum/'
#tag='continuum'
#path='/home/a.bonaldi/local2/scratch/Bonaldi/Radio_srccnt/TESTS_newtrecs/SDC2/test4/Hi/'
#tag='HI'

path='/home/a.bonaldi/local2/scratch/Bonaldi/Radio_srccnt/TESTS_newtrecs/SDC2/test4/cross/'
tag='X'
sim_side=5.4


path_dm='/data/home/a.bonaldi/local2/scratch/Bonaldi/Radio_srccnt/cones_fits/v2/'
sim_side_dm=5.

path_out='/home/a.bonaldi/local2/scratch/Bonaldi/Radio_srccnt/TESTS_newtrecs/SDC2/test4/cross_clustering/'

####end general settings


halfside_dm=sim_side_dm/2.
halfside=sim_side/2.

redshift_names=['0.01','0.02','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00','1.20','1.40','1.60','1.80','2.00','2.20','2.40','2.60','2.80','3.00','3.20','3.40','3.60','3.80','4.00','4.20','4.40','4.60','4.80','5.00','5.20','5.40','5.60','5.80','6.00','6.20','6.40','6.60','6.80','7.00','7.20','7.40','7.60','7.80','8.00','8.20','8.40','8.60','8.80','9.00','9.20','9.40','9.60','9.80','10.0']

for i in range(len(redshift_names)):


    z=redshift_names[i]
    print('********************')
    print('Processing redshift',z)
    print('********************')

    if (np.float(z) >= zmin) and (np.float(z) < zmax):
        cat_name1 = path+'catalogue_'+tag+'_z'+z+'.fits'
        cat1 = Table.read(cat_name1, format='fits')
        cat_fits1 = fits.open(cat_name1)
        #cat1 = cat_fits1[1].data
        cols1 = cat_fits1[1].columns.names

        #check if there is Mh - if not return error message

       
        #print(cols1)



        cat_name2 = path_dm+'catalogue_DM_z'+z+'.fits'
        cat2 = Table.read(cat_name2, format='fits')
        #cat_fits2 = fits.open(cat_name2)
        #cat2 = cat_fits2[1].data
        #cols2 = cat_fits2[1].columns.names

        cat_name1_out = path_out+'catalogue_'+tag+'_z'+z+'.fits' 
       
        M1 = cat1['Mh']
        print(np.min(M1),np.max(M1))

        #to deal with cross catalogues, Mh field could be flagged.  
        crosscat=0
        for k in range(len(cols1)):
            if (cols1[k] == 'Mh_1'): 
                crosscat=1 #this is a cross catalogue

        if (crosscat==1):        
            M1_bis = cat1['Mh_1']
            M1[M1==-100]=M1_bis[M1==-100.]

        print(np.min(M1),np.max(M1))

        #exit()
        M2 = cat2['Mh']
        rho2=cat2['N_2Mpc^3'] #local DM density

        #check if there is HI and if so read the MHI column
      

        ngals=len(M1)
        nhaloes=len(M2)

        attr1=np.array([M1]).T
        attr2=np.array([M2]).T

        lat1=cat1['y_coord']        #coordinates to be updated
        lon1=cat1['x_coord']        ### (default in case of no match)
        redshift1=cat1['redshift']
        
        if (crosscat==1):        
            lat1_bis = cat1['y_coord_1']
            lat1[lat1==-100]=lat1_bis[lat1==-100.]
            lon1_bis = cat1['x_coord_1']
            lon1[lon1==-100]=lon1_bis[lon1==-100.]

            redshift1_bis = cat1['redshift_1']
            redshift1[redshift1==-100]=redshift1_bis[redshift1==-100.]



        lat2=cat2['y_coord']        # halo coordinates
        lon2=cat2['x_coord']
        redshift2=cat2['redshift']

        print('number of galaxies',ngals)
        print('number of DM haloes',nhaloes)

        MHI=np.zeros(ngals)-100. #HI masses initialised as not present (-100.)

        
        #environmental dependency for HI
        HI=0 
        for k in range(len(cols1)):
            if (cols1[k] == 'MHI'): 
                HI=1

        if (HI == 1):
            print('Catalogue contains HI')
            MHI = cat1['MHI']
            rho_thr=50. #threshold on local density for HI galaxies

        #print(MHI)

        #exit()

        ###### ANALYSIS STARTS
        nsample=10 #number of nearest neighbour considered for matching
        nbrs = NearestNeighbors(n_neighbors=nsample, algorithm='kd_tree').fit(attr2)
        distances, indices = nbrs.kneighbors(attr1)

        print('number of matches',indices.shape)
        
        #print(M1[0])
        #print(M2[indices[0]])
        #print(distances[0])

        new_lat1=np.zeros(ngals) #these are the vectors to replace cat1 coordinates
        new_lon1=np.zeros(ngals)
        new_redshift1=np.zeros(ngals)


        #sort to match high mass first - more rare and difficult to match
        # go trough suggested matches and assign 
        
        indices_sortM1=np.flip(np.argsort(M1)) # indices for masses from largest to smallest
        

        for k in range(ngals):

            ###this is the galaxy we want to associate with an halo
            j=indices_sortM1[k]
            HIflag=MHI[j] # if !=-100. this is an HI source 
            #initialise new coordinates as the original ones
            new_lat1[j]=lat1[j]  
            new_lon1[j]=lon1[j]
            new_redshift1[j]=redshift1[j]
            
            #print('old coos',j,new_lat1[j],new_lon1[j])
            #proceed with clustered coordinates only in the region of overlap between cube and FoV
            if (abs(lon1[j])<=halfside_dm) and (abs(lat1[j])<=halfside_dm):
                #going in decreasing mass order, get the closest mass halo without repetition
                for jj in range(nsample):
                    
                    if (M2[indices[j,jj]] != -100.):
                        

                        new_lat1[j]=lat2[indices[j,jj]]
                        new_lon1[j]=lon2[indices[j,jj]]
                        new_redshift1[j]=redshift2[indices[j,jj]]
                        jj_match=jj  #this is to flag the mass later to avoid repetitions
                        rho=rho2[indices[j,jj]]
                        
                        break
                #in the HI case, try a different match if associated with dense environment
                #print('new coos',j,new_lat1[j],new_lon1[j])

                if (rho>rho_thr) and (HIflag != -100.): 
                    #going in decreasing mass order, get the closest mass halo with sub-threshold density and without repetition       
                    ###print('retry')
                    for jj in range(nsample):

                        if (M2[indices[j,jj]] != -100.) and (rho2[indices[j,jj]]<rho_thr):
                            new_lat1[j]=lat2[indices[j,jj]]
                            new_lon1[j]=lon2[indices[j,jj]]
                            new_redshift1[j]=redshift2[indices[j,jj]]
                            jj_match=jj  #this is to flag the mass later
                            rho=rho2[indices[j,jj]]
                            break
                #print(M1[j],M2[indices[j,jj_match]],rho)         
                M2[indices[j,jj_match]]=-100. #flagging this halo so it cannot be reused
                #print('new coos',j,new_lat1[j],new_lon1[j])
        #exit()
        print('updating coordinate information')
        cat1['y_coord']=new_lat1        # coordinates to be updated
        cat1['x_coord']=new_lon1        # (original in case of no match)
        cat1['redshift']=new_redshift1
    
        print('writing updated catalogue file')

        catout=Table()
        catout=cat1
        
        catout.write(cat_name1_out,format='fits', overwrite = True)
        
tend = time.time()
print ('...done in {0} seconds.'.format(tend-tstart))
