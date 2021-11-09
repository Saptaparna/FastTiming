import uproot3 as uproot
import awkward
import numpy as np
from scipy.sparse import csr_matrix, find
from scipy.spatial import cKDTree
from tqdm import tqdm_notebook as tqdm


import os
import os.path as osp

fname = 'NPZConversion/testNanoML_Diphoton_test.root'

print(type(fname))

test = uproot.open(fname)['Events']

import torch

from scipy.sparse import coo_matrix # to encode the cluster mappings
from sklearn.neighbors import NearestNeighbors

recHitTime = test['RecHitHGC_time'].array()
recHitEnergy = test['RecHitHGC_energy'].array()
recHitX = test['RecHitHGC_x'].array()
recHitY = test['RecHitHGC_y'].array()
recHitZ = test['RecHitHGC_z'].array()
#recHitHitR = test['RecHitHGC_hitr'].array()
recHitR = np.sqrt(recHitX*recHitX+recHitY*recHitY+recHitZ*recHitZ)
recHitTheta = np.arccos(recHitZ/recHitR)
recHitEta = -np.log(np.tan(recHitTheta/2))
print ('recHitEnergy.content = ', recHitEnergy.content)
print ('recHitR = ', recHitR)
zeros = awkward.zeros_like(recHitEta)
recHitHitR = awkward.zeros_like(recHitEta) #just setting to zero for now

features = awkward.concatenate([
            recHitEnergy,
            recHitEta,
            zeros, #indicator if it is track or not
            recHitTheta,
            recHitR,
            recHitX,
            recHitY,
            recHitZ,
            recHitTime,
            recHitHitR
            ], axis=-1) 
