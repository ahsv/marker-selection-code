# Alexander Vargo, Septermber 23 2018
# 
# Running d3e on the paul15 dataset.

##### IMPORTS
import numpy as np
import pandas as pd
import scanpy.api as sc

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80, color_map='viridis')  # low dpi (dots per inch) yields small inline figures
sc.logging.print_versions()


### USING THE PAUL15 DATA SET
adata = sc.datasets.paul15()

### READ PARAMETERS FROM FILES
# We never had this running in an array form.
ARRAY_JOB = False

import os
if ARRAY_JOB:
    # Get the array ID and switch into the correct directory for an array job
    # TODO - probably eventually use a context manager.  
    import os
    arrID = os.getenv('PBS_ARRAYID', 1)

    path = os.getcwd()
    path = path + "/" + arrID
    os.chdir(path)
    print("Current working directory: {}".format(os.getcwd()))

    # read in the max sparsity parameter from the working directory
    pFilename = "foldNum.dat"
    pFile = open(pFilename, 'r')
    foldN = pFile.readline().strip()
    foldN = int(foldN)


lookup = list(adata.obs['paul15_clusters'].cat.categories)
y = np.array([lookup.index( adata.obs['paul15_clusters'][i] ) for i in range(adata.obs['paul15_clusters'].shape[0]) ])

import shutil
from subprocess import call, Popen

# prepare the file that will be used for d3e.
# we need to replace the first line with a line consisting of the 
# groups that we want to check
def prep_file(clust,y, fname="paul-base.dat", outfname="paul-currClust.dat"):

    clustGroup = 1.0 * (y==clust)

    from_file = open(fname)
    line=from_file.readline()

    line = "GeneID"
    for group in clustGroup:
        line += "\t" + str(int(group))

    line += "\n"

    to_file = open(outfname,mode="w")
    to_file.write(line)

    shutil.copyfileobj(from_file, to_file)
    to_file.close()

    call("python D3ESplitData.py " + outfname + " 0 1 splitDir/ 150", shell=True)

#prep_file(y,0,fname="/home/ahsvargo/publicData/paul/paul-base.dat")


def D3E_for_clust(y,clust,fname="/home/ahsvargo/publicData/paul/paul-base.dat"):

    print("Working on cluster {}".format(clust), flush=True)

    # Things won't work for an array job right now.
    if ARRAY_JOB:
        dataname = "/home/ahsvargo/publicData/paul/fold" + str(foldN) + "/paul-currClust.dat"
        outname = "/home/ahsvargo/publicData/paul/fold" + str(foldN) + "/paul-result.out"
    else:
        dataname = "/home/ahsvargo/publicData/paul/paul-currClust.dat"
        outname = "/home/ahsvargo/publicData/paul/paul-result"
        splitpath = "/home/ahsvargo/publicData/paul/D3E_tests/splitDir/"

    prep_file(clust,y,fname=fname, outfname=dataname)

    files = os.listdir(splitpath)
    plist = []
    for name in files:
        realName, ext = os.path.splitext(name)
        outfile = outname + realName  + ".out"
        print(outfile)

        p = Popen("python /home/ahsvargo/publicData/paul/D3E_tests/D3ECmd.py " + 
                splitpath + name + " " + outname + 
                " 1 0 -m 0 -t 0 -z 0 -n 1 -v", shell=True)

        plist.append(p)

    exit_codes = [p.wait() for p in plist]
    return exit_codes

D3E_for_clust(y,0)

