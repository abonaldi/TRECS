from sklearn.neighbors import NearestNeighbors
import numpy as np
import numpy as np
import os, glob, sys
from astropy.io import fits
import matplotlib
from matplotlib import pyplot as plt
import astropy
from collections import Counter
from astropy.cosmology import LambdaCDM
import time
from astropy.table import Table
from astropy.table import hstack
from astropy.table import vstack

tstart = time.time()
compare_cats = 0
line_flux_cut = 1.6 #Jy Hz


def match(cat1, cat2,path,write, properties=None):
    if not os.path.isfile(cat1):
        cat1 = cat1.replace(cat1.split('_')[-1], 'z0.20.fits') # use as dummy to get columns
        HI = False
       
       
    else:
        HI = True    

 
   
        
          
    cat_name_HI = cat1
    cat_name_cont = cat2

    cat_fits_HI = fits.open(cat_name_HI)
    cat_fits_cont = fits.open(cat_name_cont)

    cat_HI = cat_fits_HI[1].data
    cat_cont = cat_fits_cont[1].data

    cols_HI = cat_fits_HI[1].columns.names
    cols_cont = cat_fits_cont[1].columns.names

    cat_HI_table = Table.read(cat_name_HI, format='fits')
    cat_cont_table = Table.read(cat_name_cont)


    


    for i in range(len(cols_cont)):
        if cols_cont[i] in cols_HI:
            cols_cont[i] = cols_cont[i]+'_1'
  


    

    # how to convert a recarray or fits table to np array:
    cat_HI_np = np.array(cat_HI_table).view(np.float32).reshape((np.array(cat_HI_table).shape + (-1,)))
    cat_cont_np = np.array(cat_cont_table).view(np.float32).reshape((np.array(cat_cont_table).shape + (-1,)))



    print (cols_cont, cols_HI)

    if HI:
        print ('cat lengths', cat1.split('/')[-1],  len(cat_HI), len(cat_cont))

        MHI_HI = cat_HI['MHI']
        MH_HI = cat_HI['Mh']
        #r_HI_HI = cat_HI['HI size']
        line_flux_HI = cat_HI['HI flux']/1000 # from mJy to Jy
        incl_HI = cat_HI['inclination']
        z_HI = cat_HI['redshift']
        OptClass_HI = cat_HI['OptClass']
    

        MHI_cont = cat_cont['MHI_pred']
        MH_cont = cat_cont['Mh_1']
        #r_HI_cont = cat_cont['HI size']
        incl_cont = cat_cont['inclination_1']
        z_cont = cat_cont['redshift_1']
        mass_function = 0


        if mass_function:   
            cont_optclasses =[ MHI_cont[cat_cont['RadioClass']==1],MHI_cont[cat_cont['RadioClass']>3]]
            labels = ['late-type', 'late-type + AGN']
            colors = ['red', 'pink']

            plt.clf()

            norm = False
            #plt.hist(MHI_HI,range = (7.5, 12), bins = 100,log=True, histtype='step', fill=False,label = 'MHI', alpha = 1, normed=norm)
            #plt.hist(cont_optclasses,range = (7.5, 12),stacked = True, histtype='step', fill=False,bins = 100, log=True, alpha = 1, normed=norm, color = colors, label = labels)
          


            plt.legend()
            plt.xlabel(r'log MHI (M$_{\odot}$)')
            plt.ylabel('N sources')
            plt.title(cat2.split('_')[-1].split('.fits')[0])
            plt.savefig('cross/plots/HI_number_counts%s.png'%cat2.split('_')[-1].split('.fits')[0])
            return




        #work out line_flux_pred from MHI_pred and dont match any continuum sources with line_flux_pred below line flux count


        H=67.0
        M=0.32
        L=0.68
        c = 2.99792458e8
        G = 6.67408e-11
        cosmo = LambdaCDM(H0 = H, Om0 = M, Ode0 = L)
        D_L_cont = cosmo.luminosity_distance(z_cont).value # Mpc
        
        line_flux_cont = 10**MHI_cont/ (49.8 * D_L_cont**2)





        print (MHI_cont)
        print (len(cat_cont), 'continuum sources')
        print (len(cat_cont[line_flux_cont>=line_flux_cut]), 'continuum sources predict HI flux above HI cut')
        print (len(cat_cont[line_flux_cont<line_flux_cut]), 'continuum sources will not be matched with HI')
        print (len(cat_HI), 'HI sources')
        print (len(cat_HI)-len(cat_cont[line_flux_cont>=line_flux_cut]), 'HI sources will not be matched with continuum')
        print (len(cat_cont)+len(cat_HI)-len(cat_cont[line_flux_cont>=line_flux_cut]), 'unique sources in catalogue')


        unmatched_cont = cat_cont_np[line_flux_cont<line_flux_cut]

        unmatched_cont_empty = np.zeros((unmatched_cont.shape[0],cat_HI_np.shape[1]))-100
 
        unmatched_cont_stack = np.hstack((unmatched_cont_empty, unmatched_cont))
   
        matched_cat_cont_np = cat_cont_np[line_flux_cont>=line_flux_cut]

        # find lowest N MHI sources in HI cat, where N is the number of surplus HI sources after matching
        # with all continuum sources with predicted flux over HI flux threshold

        N_unmatched_HI = len(cat_HI)-len(cat_cont[line_flux_cont>=line_flux_cut]) 
        print ('N_unmatched_HI',N_unmatched_HI)
        # catch values less than zero
        N_unmatched_HI=np.max((N_unmatched_HI, 0))
        print ('N_unmatched_HI',N_unmatched_HI)

        print (line_flux_cont.shape)
        print (line_flux_HI.shape)
        # value of MHI_HI of Nth lowest source after sorting in order of MHI_HI
        sorted_line_flux_HI = np.sort(line_flux_HI)
        HI_cat_line_flux_cut = sorted_line_flux_HI[N_unmatched_HI]
        print ('all HI sources with line flux below', HI_cat_line_flux_cut, 'Jy will not be matched')





        unmatched_HI = cat_HI_np[line_flux_HI<HI_cat_line_flux_cut]

        unmatched_HI_empty = np.zeros((unmatched_HI.shape[0], cat_cont_np.shape[1]))-100

        unmatched_HI_stack = np.hstack((unmatched_HI, unmatched_HI_empty))

        matched_cat_HI_np = cat_HI_np[line_flux_HI>=HI_cat_line_flux_cut]


        all_cols = cols_HI+cols_cont

        unmatched_HI_table = Table()
        for i, col in enumerate(all_cols):
            unmatched_HI_table[col] = unmatched_HI_stack[:,i]


        unmatched_cont_table = Table()
        for i, col in enumerate(all_cols):
            unmatched_cont_table[col] = unmatched_cont_stack[:,i]            




        matched_MHI_HI = MHI_HI[line_flux_HI>=HI_cat_line_flux_cut]
        matched_MHI_cont = MHI_cont[line_flux_cont>=line_flux_cut]
        print (matched_MHI_HI.shape, matched_MHI_cont.shape)


        sorted_matched_MHI_HI = matched_MHI_HI[np.argsort(matched_MHI_HI)]

        sorted_matched_MHI_cont = matched_MHI_cont[np.argsort(matched_MHI_cont)]
        for i in range(len(sorted_matched_MHI_HI)):
            print (sorted_matched_MHI_HI[i], sorted_matched_MHI_cont[i])


        print (sorted_matched_MHI_cont.shape)
        print(sorted_matched_MHI_HI.shape)
        both = np.vstack((sorted_matched_MHI_HI, sorted_matched_MHI_cont))
        print (both)

        matched_cat_HI_np = matched_cat_HI_np[np.argsort(matched_MHI_HI)]
        newcat2 = matched_cat_cont_np[np.argsort(matched_MHI_cont)]
        




        # now only need to match the matched catalogues, and reserve the unmatched to stack at the end

   


        matching_cat_HI_table = Table()
        for i, col in enumerate(cols_HI):
            matching_cat_HI_table[col] = matched_cat_HI_np[:,i]

    elif not HI:
                # make a numpy array for now
       
        newcat1 = np.zeros((cat_cont_np.shape[0], cat_HI_np.shape[1]))-100   
        matching_cat_HI_table = Table()
        for i, col in enumerate(cols_HI):
            matching_cat_HI_table[col] = newcat1[:,i]


        newcat2 = cat_cont_np

    # might need to make HI table from np array here as it is reordered

    # make it into a fits table
    cat = Table()
    for i, col in enumerate(cols_cont):
        cat[col] = newcat2[:,i]

    t_new = hstack([matching_cat_HI_table, cat])


    plt.clf()

    plt.scatter(t_new[t_new['OptClass_1']==2]['MHI'], t_new[t_new['OptClass_1']==2]['MHI_pred'], label = 'spirals')
    plt.scatter(t_new[t_new['OptClass_1']==1]['MHI'], t_new[t_new['OptClass_1']==1]['MHI_pred'], label = 'ellipticals')
    plt.xlabel(r'log MHI (M$_{\odot}$)')
    plt.ylabel(r'log MHI_pred (M$_{\odot}$)')
    plt.legend()
    plt.savefig(path+'plots/MHI_pred_vs_MHI%s.png'%cat2.split('_')[-1].split('.fits')[0])


    plt.clf()

    plt.scatter(t_new[t_new['OptClass_1']==2]['MHI'], t_new[t_new['OptClass_1']==2]['MHI']-t_new[t_new['OptClass_1']==2]['MHI_pred'], label = 'spirals')
    plt.scatter(t_new[t_new['OptClass_1']==1]['MHI'], t_new[t_new['OptClass_1']==1]['MHI']-t_new[t_new['OptClass_1']==1]['MHI_pred'], label = 'ellipticals')
    plt.xlabel(r'log MHI (M$_{\odot}$)')
    plt.ylabel(r'log MHI - log MHI_pred (M$_{\odot}$)')
    plt.legend()
    plt.savefig(path+'plots/MHI_pred_vs_MHI_res%s.png'%cat2.split('_')[-1].split('.fits')[0])

  
    if HI:

        # vstack the unmatched here
        t_new_all = vstack([t_new, unmatched_HI_table, unmatched_cont_table])    
    else:
        t_new_all = t_new


    plot_MHI_dist = 0
    if plot_MHI_dist:
        plt.clf()
        plt.hist(t_new[t_new['OptClass_1']==2]['MHI_pred'], alpha = 0.5, label = 'spirals')
        plt.hist(t_new[t_new['OptClass_1']==1]['MHI_pred'], alpha = 0.5, label = 'ellipticals')
        plt.xlabel(r'log MHI_pred (M$_{\odot}$)')
        plt.ylabel('N')
        plt.legend()
        plt.savefig('MHI_pred%s.png'%cat2.split('_')[-1].split('.fits')[0])
        plt.clf()

       

        plt.hist(t_new[(t_new['RadioClass']==1) ]['MHI_pred'], alpha = 0.5, label = 'SFG-late')
        plt.hist(t_new[(t_new['RadioClass']==2 ) ]['MHI_pred'], alpha = 0.5, label = 'SFG-early')
        plt.hist(t_new[(t_new['RadioClass']>3 )]['MHI_pred'], alpha = 0.5, label = 'AGN')   
        plt.xlabel(r'log MHI_pred (M$_{\odot}$)')
        plt.ylabel('N')
        plt.legend()
        plt.savefig('MHI_pred_radioclass%s.png'%cat2.split('_')[-1].split('.fits')[0])
        plt.clf()




    outf = (path+cat2.split('/')[-1].replace("continuum", "X"))
    print ('writing to.. ', outf)
    if write:
        t_new_all.write(outf,format='fits', overwrite = True)



if __name__=='__main__':

    path_c='/home/a.bonaldi/SDC2/TRECS_outputs/sdc2_fullcube2/continuum/'

    path_hi='/home/a.bonaldi/SDC2/TRECS_outputs/sdc2_fullcube2/HI/'

    path_cross='/home/a.bonaldi/SDC2/TRECS_outputs/sdc2_fullcube2/cross/' 


    #redshift_names=['0.01','0.02','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00','1.20','1.40','1.60','1.80','2.00','2.20','2.40','2.60','2.80','3.00','3.20','3.40','3.60','3.80','4.00','4.20','4.40','4.60','4.80','5.00','5.20','5.40','5.60','5.80','6.00','6.20','6.40','6.60','6.80','7.00','7.20','7.40','7.60','7.80','8.00','8.20','8.40','8.60','8.80','9.00','9.20','9.40','9.60','9.80','10.0']

    cat_name_pattern =  path_c+'catalogue_continuum*.fits'


    continuum_files = glob.glob(cat_name_pattern)
    
    print (continuum_files)

    for i in continuum_files:
        cat_name_cont = i
      
        cat_name_HI =  i.replace('continuum', 'HI') # '/Users/p.hartley/Dropbox (SKA)/science_team/HI_sims/code/SDC2/catalogues/hi/catalogue_HI_z0.%02d.fits'%i
        print (cat_name_cont)
    
        match(cat_name_HI,cat_name_cont, path_cross,write = True)
        


    tend = time.time()
    print ('...done in {0} seconds.'.format(tend-tstart))
