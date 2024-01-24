#@PYTHONSHEBANG@
# 
# overwrites the sky coordinates (x_coord, y_coord, redshift) of galaxies with those of a lightcone of a cosmological simulation, thus reproducing clustering
# matching galaxies between the DM catalogue and the observed catalogue based on DM halo mass
# use nearest neighbours for associating mass
# if HI is present, use local density as well to introduce environmental dependencies (HI tends to avoid high local density)
# 
# to improve speed, only around 1000 haloes per galaxy in the same mass bin are retained in the DM catalogue.  
#####


# History
#---------------------------------
# A. Bonaldi 14/7/21 first version
# T. Ronconi 22/3/23 I/O mods and argument parsing





from sklearn.neighbors import NearestNeighbors
import numpy as np
import os, glob, sys
from astropy.io import fits
import astropy
from astropy.table import Table

import time
 
tstart = time.time()

################################################################
# Read parameter file

# Check a parameter file is passed
if len(sys.argv) != 2 :
    raise RuntimeError(f'script {sys.argv[0]} needs a parameter file to run')

# check that the path passed is a file
# and in case read it
if os.path.isfile(sys.argv[1]) :
    with open(sys.argv[1], 'r') as f :
        lines = f.readlines()
else :
    raise RuntimeError(f'provided path {sys.argv[1]} is not a file')

# remove commented lines and \newlines
lines = [ line[:-1].split('#')[0] for line in lines if line[0] != '#' ]

# create content dictionary from parameter file
content = { k.lstrip().rstrip() : v.lstrip().rstrip()
            for (k,v) in [
                    line.split("=")
                    for line in lines
                    if len(line) > 0
            ]
}

################################################################
# Parse arguments

zmin = float(content['z_min'])
zmax = float(content['z_max'])
print( f'Cross matching in the redshift range {zmin:.3f} <= z <= {zmax:.3f}' )

# Field of view of the crossmatched catalogue
# this should be <= than the catalogue without clustering.
# If > than the DM catalogue, clustering is done only in the central part
fov = float(content['sim_side'])

# try reading tag from parameter file
# (if not provided set to default)
tag = ''
try :
    tag = content['tag_DMmatch']
except KeyError :
    tag = 'continuum'

# output directory parsing and check existence
if not os.path.isdir(content['outdir']) :
    raise RuntimeError(f'output directory {content["outdir"]} does not exist')
path_out=os.path.join( content['outdir'], '_'.join( [ 'raw', tag, 'clustered' ] ) )
if not os.path.isdir(path_out) :
    try :
        os.mkdir(path_out)
    except :
        raise

# input directory containing data to cross-match
path = os.path.join(content['outdir'], 'raw_'+tag)
if not os.path.isdir(path) :
    raise RuntimeError(f'input continuum directory {path} does not exist')

# input directory with the clustering properties
path_dm = os.path.join( content['TRECSinputs'], 'TRECS_Inputs', 'DMhalo_catalogues' )
if not os.path.isdir(path_dm) :
    raise RuntimeError(
        f'input continuum directory {path_dm} does not exist, '
        'did you download the inputs database?'
    )

# shifts from centre of simulation to centre of field
x_shift = None
try :
    x_shift = float( content['x_shift'] )
except KeyError :
    x_shift = 0.0
except :
    raise
y_shift = None
try :
    y_shift = float( content['y_shift'] )
except KeyError :
    y_shift = 0.0
except :
    raise

# shifts from centre of simulation to centre of field 
x_shift_dm = None
try :
    x_shift_dm = float( content['x_shift_dm'] )
except KeyError :
    x_shift_dm = 0.0
except :
    raise
y_shift_dm = None
try :
    y_shift_dm = float( content['y_shift_dm'] )
except KeyError :
    y_shift_dm = 0.0
except :
    raise

# Read the general random seed and initialize random number generator:
seed = None
try :
    seed = 888 * int( content['seed'] )
except KeyError :
    pass
except :
    raise

#rngen = np.random.default_rng( seed = seed )
np.random.seed( seed = seed )

################################################################
# Hard-coded settings

# size of square simulation (degs)
sim_side=fov

# size of square DM simulation (degs)
sim_side_dm=5.

# the number of objects above which the
# thinning speed-up technique is applied
# when cross-matching
Ncrit_thinning = 10000

########################################### end general settings

halfside_dm=sim_side_dm/2.
halfside=sim_side/2.

# if we simulated a smaller FoV than the dm cone, we select just a portion of the cone
if (halfside_dm > fov/2.):
    halfside_dm=fov/2.

if (halfside > fov/2.):
    halfside=fov/2.

#check that the shift does not go outside the catalogue areas
if (halfside_dm+np.abs(x_shift_dm) > sim_side_dm/2.):
    print('Shift on x for DM cat too big')
    exit()

if (halfside_dm+np.abs(y_shift_dm) > sim_side_dm/2.):
    print('Shift on y for DM cat too big')
    exit()

if (halfside+np.abs(x_shift) > sim_side/2.):
    print('Shift on x for cat too big')
    exit()

if (halfside+np.abs(y_shift) > sim_side/2.):
    print('Shift on y for cat too big')
    exit()


redshift_names=['0.01','0.02','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00','1.20','1.40','1.60','1.80','2.00','2.20','2.40','2.60','2.80','3.00','3.20','3.40','3.60','3.80','4.00','4.20','4.40','4.60','4.80','5.00','5.20','5.40','5.60','5.80','6.00','6.20','6.40','6.60','6.80','7.00','7.20','7.40','7.60','7.80','8.00','8.20','8.40','8.60','8.80','9.00','9.20','9.40','9.60','9.80','10.0']

mask_z = np.zeros_like(redshift_names, dtype=bool)
for i in range(len(redshift_names)):
    
    z=redshift_names[i]

    if (np.float(z) >= zmin) and (np.float(z) <= zmax):
        print('********************')
        print('Processing redshift',z)
        print('********************')
        cat_name1 = os.path.join( path, 'catalogue_'+tag+'_z'+z+'.fits' )
        # nothing done if catalogue does not exist
        if (os.path.isfile(cat_name1) == True):

            cat1 = Table.read(cat_name1, format='fits')
            cat_fits1 = fits.open(cat_name1)
            cols1 = cat_fits1[1].columns.names


            #select FoV only
            cat1=cat1[(np.abs(cat1['x_coord']+x_shift) <= halfside)*(np.abs(cat1['y_coord']+y_shift) <= halfside)]     

        
            cat_name2 = os.path.join( path_dm, 'catalogue_DM_z'+z+'.fits' )
            cat2 = Table.read(cat_name2, format='fits')

            #here select only the portion of the catalogue that we need
            cat2=cat2[(np.abs(cat2['x_coord']+x_shift_dm) <= halfside_dm)*(np.abs(cat2['y_coord']+y_shift_dm) <= halfside_dm)] 

       
            cat_name1_out = os.path.join( path_out, 'catalogue_'+tag+'_clustered_z'+z+'.fits' )
       

            #start reading catalogue values
            M1 = cat1['Mh']
            M2 = cat2['Mh']
        

            #check if there is HI and if so read the MHI column
            
            #here I could take a portion of the DM catalogue for low masses to speed-up things 

            #apply a mass cut to DM catalogue 
            minmass=np.min(M1) #this is the minimum mass to be matched, decreased to be conservative
            cat2=cat2[cat2['Mh'] >= minmass]
            M2 = cat2['Mh']
            rho2=cat2['N_2Mpc^3'] #local DM density

            #apply thinning of DM catalogue to speed-up the process
            if M2.size > Ncrit_thinning or M1.size > Ncrit_thinning :
            
                #histogram of M1 and M2 to see and work out ratio between galaxies and eligible haloes
                maxmass=np.max((np.max(M1),np.max(M2)))
                nbins=50
                d_m1=np.histogram(M1,range=(minmass,maxmass),bins=nbins)
                d_m2=np.histogram(M2,range=(minmass,maxmass),bins=nbins)

                #print('dm1',d_m1[0])
                #print('dm2',d_m2[0])

                #       TODO: avoid the case where d_m1=0
                ratio=(d_m2[0]/(d_m1[0]+1)) #this is how many dark haloes per observable galaxy as a finction of mass

                binning=d_m2[1]
                
                
                thr=np.zeros(len(ratio))+1.
                thr[(ratio>1.)]=1./ratio[ratio>1.]*10000.

                select=np.zeros(len(M2), dtype=bool) #vector to record the portion of the catalogue to keep for crossmatching
                
                
                rnd=np.random.uniform(low=0, high=1, size=len(M2))
                
       
                for k in range(len(M2)):
                    value=M2[k]
                    idx = (np.abs(binning - value)).argmin()
                    if (idx>=nbins):
                        idx=nbins-1

                    if (rnd[k]<=thr[idx]):
                        select[k]=True

            
                cat2=cat2[select] #selection of eligible haloes complete
                M2 = M2[select]
                rho2 = rho2[select]
#                d_m2_new=np.histogram(M2,range=(minmass,maxmass),bins=nbins)

#                print(d_m2_new[0])
#                print(np.sum(d_m2_new[0]))
                
                print('Number of original DM haloes',

len(rnd))
                print('Number of sub-sampled DM haloes (the 2 should be the same)',
                      np.sum(select),len(cat2))

                print('Number of galaxies',
                      len(cat1))



            ######## end thinning for speed-up
    
            # M2 = cat2['Mh']
            # rho2=cat2['N_2Mpc^3'] #local DM density
        
            ngals=len(M1)
            nhaloes=len(M2)
            if ngals == 0 :
                print( "No galaxies in current redshift slice: skipping" )
                break

            if nhaloes != 0:

                attr1=np.array([M1]).T   #observed catalogue
                attr2=np.array([M2]).T   #lightcone
                
                lat1=cat1['y_coord']+y_shift        #coordinates to be updated
                lon1=cat1['x_coord']+x_shift        ### (default in case of no match)
                redshift1=cat1['redshift']
        
                print('Cat1 selected coordinates')
                print(np.min(lon1),np.max(lon1))
                print(np.min(lat1),np.max(lat1))
        

                lat2=cat2['y_coord']+y_shift_dm        # halo coordinates
                lon2=cat2['x_coord']+x_shift_dm
                
                print('Cat2 selected coordinates')
                print(np.min(lat2),np.max(lat2))
                print(np.min(lon2),np.max(lon2))
        
                redshift2=cat2['redshift']

        
                print('number of galaxies',ngals)
                print('number of DM haloes',nhaloes)

                MHI=np.zeros(ngals)-100. #HI masses initialised as not present (-100.)

        
                #environmental dependency for HI
                #check if HI is present in the column header
                HI=0
                rho_thr=50. #threshold on local density for HI galaxies
                for k in range(len(cols1)):
                    if (cols1[k] == 'MHI'): 
                        HI=1

                if (HI == 1):
                    print('Catalogue contains HI field')
                    MHI = cat1['MHI']
                


                ###### ANALYSIS STARTS
                print('Start nearest neighbour')
        
                nsample=np.min([20,nhaloes])
                nbrs = NearestNeighbors(n_neighbors=nsample, algorithm='kd_tree').fit(attr2)
                distances, indices = nbrs.kneighbors(attr1)
                print('Neaerst neighbour done')
        
                print('number of matches',indices.shape)
        
                new_lat1=np.zeros(ngals) #these are the vectors to replace cat1 coordinates
                new_lon1=np.zeros(ngals)
                new_redshift1=np.zeros(ngals)


                #sort to match high mass first - more rare and difficult to match
                # go trough suggested matches and assign 
        
                indices_sortM1=np.flip(np.argsort(M1)) # indices for observed masses from largest to smallest
        

                for k in range(ngals):

                    ###this is the galaxy we want to associate with an halo
                    j=indices_sortM1[k]
                    HIflag=MHI[j] # if !=-100. this is an HI source 

                    #initialise new coordinates as the original ones
                    new_lat1[j]=lat1[j]  
                    new_lon1[j]=lon1[j]
                    new_redshift1[j]=redshift1[j]
            
                
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
                    

                        if (rho>rho_thr) and (HIflag != -100.): #if the galaxy has HI conterpart and the region is dense:
                            # try match again with additional constraint on local density
                            for jj in range(nsample):

                                if (M2[indices[j,jj]] != -100.) and (rho2[indices[j,jj]]<rho_thr):
                                    new_lat1[j]=lat2[indices[j,jj]]
                                    new_lon1[j]=lon2[indices[j,jj]]
                                    new_redshift1[j]=redshift2[indices[j,jj]]
                                    jj_match=jj  #this is to flag the mass later
                                    rho=rho2[indices[j,jj]]
                                    break
      
                        M2[indices[j,jj_match]]=-100. #flagging this halo so it cannot be reused
                        
                        ### end for k in range(ngals):
      
                # update coordinates and redshift after associating with DM haloes
                print('updating coordinate information')
                cat1['y_coord']=new_lat1        # coordinates to be updated
                cat1['x_coord']=new_lon1        # (original in case of no match)
                cat1['redshift']=new_redshift1

                ############ end of nearest neighbours
            
            print('writing updated catalogue file')

            # for final slices file
            mask_z[i] = True

            catout=Table()
            catout=cat1
        
            catout.write(cat_name1_out,format='fits', overwrite = True)


# Writing redshift slices file
with open(
        os.path.join(
            content['outdir'],
            '_'.join(['slices',
                      tag,
                      'clustered'])+'.dat'
        ), 'w' ) as outf :
    for z in np.array( redshift_names )[mask_z] : outf.write(f'{z}\n')    
        
tend = time.time()
print ('...done in {0} seconds.'.format(tend-tstart))
