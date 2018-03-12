r"""
Module: combine.py

the module combines pickle files from a kMC
simulation.

See Also
--------
kmc.save_pickle() 

"""

import numpy as np
import os
import os.path as path
import types
import pickle

def combine_pickles(outfilename='transitions_combined.pickle'):
    r"""Combines all pickle files in current directory.
    
        The method identifies all pickle files in directory,
        prepares lists and dictionaries to combine the pickles,
        and finally writes a file with filename 'outfilename'
        to disc.
        
        Parameters
        ----------
        outfilename : str
            path to output the combined pickle file.
            
     """
    # Find all picles in current dir
    files_tmp = []
    for file in os.listdir('.'):
        if file.endswith('.pickle') and 'combined' not in file \
            and outfilename != file:
            
            files_tmp.append(file)

    # Sort pickle files by number:
    Numbers = []
    for i in range(len(files_tmp)):
        Numbers.append(int(files_tmp[i].split('step')[1].\
                       split('.')[0]))
    
    srt = np.argsort(Numbers)

    files = [files_tmp[srt[j]] for j in range(len(Numbers))]
 
    d1 = pickle.load(open(files[0],'r'))
    keylist = d1.keys()
    out = {}
    # Prepare out dict to be filled.
    for k in keylist:
        out[k] = {}
        if k == 'Siteids' or  k== 'stypes' or  k=='Nsites' or\
           k=='atoms' or k =='parameters':
            out[k] = d1[k]

     
    # Prepare stype_ev and other
    out['stype_ev'] = {}
    for u in d1['stype_ev'].keys():
        out['stype_ev'][u] = np.zeros(len(d1['stype_ev'][u]))

    out['stype_ev_other'] = {}
    for u in d1['stype_ev_other'].keys():
        out['stype_ev_other'][u] = np.zeros(len(\
                                    d1['stype_ev_other'][u])) 

    # Prepare lists:
    out['time'] = []
    out['covered'] = []
    out['nevents'] = np.zeros(len(d1['nevents']))

    # Now load each file and run through its keys.
    for i,f in enumerate(files):
        dcur = pickle.load(open(f,'r'))
        for k in keylist:
        
            if k == 'stype_ev':  
                for v in dcur['stype_ev'].keys():
                    out['stype_ev'][v] += np.array(
                                        dcur['stype_ev'][v])

            elif k == 'stype_ev_other':
                for v in dcur['stype_ev_other'].keys():
                    out['stype_ev_other'][v] += np.array(
                                  dcur['stype_ev_other'][v])

            elif k == 'time': 
                out['time'].extend(dcur['time'])

            elif k == 'covered':
                out['covered'].extend(dcur['covered'])
     

            elif k == 'nevents':   
                out['nevents'] += np.array(dcur['nevents'])


    f = open(outfilename,'w')
    pickle.dump(out,f)
    f.close()


def combine_files_iso():
    r"""Combines pickles in all dirs named run_* in turn.
    
        Method cd's into all directories named 'run_X', and
        calls combine_pickles. The result in each directory
        is stored in dir: ./combined as 'rX'.
        
        Method is relevant when multiple parallel simulations
        have been run in parallel, and error estimates are to
        be made between them.
    
     """

    # List the directories:
    fld_tmp = []
    for ent in os.listdir('.'):
        if path.isdir(ent) and 'run' in ent :
            fld_tmp.append(ent)
    
    # Make combined dir.
    if not path.isdir('combined'):
        os.mkdir('combined')
    # Cd into each and combine their individual pickles
    for f in fld_tmp:
	print 'Combining ', f
        os.chdir(f)        
        fldnumber = f.split('_')[1]
        if not path.exists('../combined/r'+fldnumber+'.pickle'):
        	combine_pickles(outfilename='../combined/r'+fldnumber+'.pickle')
        
        os.chdir('..')

if __name__=='__main__':
    combine_files_iso()
