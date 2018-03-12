r"""
Module: analysis.py

module contains methods to analyze data
generated from MonteCoffee kMC runs.

See Also
--------
kmc.save_pickle() 

"""

import pickle
import os
import numpy as np
from combine import combine_pickles
from user_sites import Site



# Fraction of points to include in analysis:
fracinc = 1/3. 

def event_histogram():
    r"""Prints the rates of different events on the 
        different site types.
        
        Method loads in all pickle files in directory './', 
        and calculates averages rates bewtween all pickle 
        files in './'.
    
        See Also
        --------
        combine.py
        kmc.save_pickle()
     
    """
    
    print 'Analyzing Events Frequencies'
    
    # Load each pickle file and add to the lists.
    files_tmp = []
    for file in os.listdir('.'):
        if file.endswith('.pickle'):
            files_tmp.append(file)
    
    # Prepare the data lists:
    with open(files_tmp[0],'r') as f:
        d1 = pickle.load(f)
        stypes = d1['stypes']   
        stype_COads = {}
        stype_COads_other = {}
        stype_COdes = {}
        stype_COdes_other = {}
        for k in list(set(stypes)):
            stype_COads[k] = []
            stype_COads_other[k] = []
            stype_COdes[k] = []
            stype_COdes_other[k] = []

     

    for j, pickle_file in enumerate(files_tmp):
        with open(pickle_file,'r') as f:
            print 'reading file ', pickle_file
            d = pickle.load(f)
            time = d['time']
            stype_ev = d['stype_ev']
            stype_ev_other = d['stype_ev_other']
           
            for k in stype_ev.keys():
                stype_COads[k].append((stype_ev[k][0]+
                               stype_ev_other[k][0])/time[-1]) 

                stype_COdes[k].append((stype_ev[k][1]+
                               stype_ev_other[k][1])/time[-1])

                
    
         
    # Take averages of all .pickle files:
    for k in stype_COads.keys():
        print 'Stype ', k,'\n', '-'*20, '\n' 
        print 'CO adsorption rate : ',\
               np.mean(stype_COads[k])
        print 'CO desorption rate : ',\
               np.mean(stype_COads[k]), '\n', '-'*20,'\n'
   


def analyze_coverages():
    r"""Analyzes the average coverage in a kMC run.
    
        Method should be called in a directory with
        pickle files combined using the combine.py
        module. 

        Occupations are read into lists together with
        the time and site-types. 

        See Also
        --------
        combine.py
        kmc.save_pickle()
     
    """

    # Load each pickle file and add to the lists.
    covs_CO = []
    sd_CO = []
    files_tmp = []
    for file in os.listdir('.'):
        if file.endswith('.pickle') and 'transitions' in file:
            files_tmp.append(file)
    
    for fcur in files_tmp:
         # Load pickle
        with open(fcur,'r') as f:
            d = pickle.load(f)
    
        time = d['time']
        stypes = d['stypes']
        covs = d['covered']
        mcstep = d['mcstep']

        # Include the last 2/3 of the simulation. 
        points = range(int(len(covs)*fracinc),len(covs))
        cov_stype_CO = [np.zeros(max(stypes)+1) for n in points]
    
        # Find the number of sites with each site-type(stype)
        Nstypes = np.zeros(max(stypes)+1)
        for i in range(len(stypes)):
            Nstypes[stypes[i]] += 1.

        # For each point in time:
        for step in range(len(points)):
            # For each site type
            for i in range(len(stypes)): 
                if covs[step][i] == 1:
                    cov_stype_CO[step][stypes[i]] += 1./Nstypes[stypes[i]]
            
        # CO on stype 0:
        cov_CO = [cov_stype_CO[i][1] for i in range(len(cov_stype_CO))]
        covs_CO.append(np.mean(cov_CO[0])) 
        sd_CO.append(np.std(cov_CO[0]))

     

    print 'Average CO coverages Cov, STD\n'+'-'*20+'\n'
    print 'CO cov :', np.mean(covs_CO), np.std(covs_CO)
    print '\n'+'-'*20+'\n'



def coverage_single(pickle_file):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.rcParams['axes.linewidth'] = 3 #set the value globally
    mpl.rcParams['mathtext.default']='rm'
    mpl.rcParams['mathtext.rm'] = 'Arial'
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rc("font",size=22)
    sz=22
    lw=2.5
    with open(pickle_file,"r") as f:
        d = pickle.load(f)
    covs = d['covered']
    time = d['time']
    stypes = d['stypes']

    cov_111 = np.zeros(len(covs))
    cov_100 = np.zeros(len(covs))
    cov_edge = np.zeros(len(covs))
    cov_cnr = np.zeros(len(covs))
    N111 = float(len([s for s in stypes if s == 3]))
    N100 = float(len([s for s in stypes if s == 2]))
    Nedge = float(len([s for s in stypes if s == 1]))
    Ncnr = float(len([s for s in stypes if s == 0]))


    for step in range(len(covs)):
        for sid,st in enumerate(stypes):
            if covs[step][sid] == 1:            
                if st == 3:
                    cov_111[step] += 1./N111
                elif st == 2:
                    cov_100[step] += 1./N100
                elif st == 1:
                    cov_edge[step] += 1./Nedge
                elif st == 0:
                    cov_cnr[step] += 1./Ncnr


    colbulk = np.array([256,256,256])/256.
    col100 = np.array([178,171,210])/256.
    col111 = np.array([94,60,153])/256.
    coledge = np.array([253,184,99])/256.
    colcnr = np.array([230,97,1])/256.


    plt.xlabel(r"Time (10$^{-6}$ s)",fontsize=sz)
    plt.ylabel("Site Coverage",fontsize=sz)
    plt.ylim((-0.05,1.05))
    plt.tick_params(axis='both', labelsize=sz,length=10, width=2, which='major')
    
    plt.plot(np.array(time)*1E6,cov_100,'-',lw=lw,color=col100)
    plt.plot(np.array(time)*1E6,cov_edge,'-',lw=lw,color=coledge)
    plt.plot(np.array(time)*1E6,cov_cnr,'-',lw=lw,color=colcnr) 
    plt.plot(np.array(time)*1E6,cov_111,'-',lw=lw,color=col111)   
    plt.savefig("covs.svg")
    plt.show()

if __name__ == '__main__':

    # Average over multiple runs :
    #analyze_coverages()
    #event_histogram()
    coverage_single("r0.pickle")
