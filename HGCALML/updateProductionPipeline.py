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
recHitHitR = awkward.zeros_like(recHitEta) #just setting to zero for now and acknowledging existence of branch

recHitFeatures = awkward.concatenate([
            recHitEnergy,
            recHitEta,
            zeros, #indicator if it is track or not
            recHitTheta,
            recHitR,
            recHitX,
            recHitY,
            recHitZ,
            recHitTime
            ], axis=-1)

##taken from Thomas' code here: https://github.com/cms-pepr/pytorch_cmspepr/blob/main/torch_cmspepr/dataset.py#L104-L170 with minor updates

recHitTruthEnergy = test['MergedSimCluster_boundaryEnergy'].array()


class SinglePhotonDataset(Dataset):
    """SinglePhoton dataset.

    Features in x:
    0 recHitEnergy,
    1 recHitEta,
    2 zeroFeature, #indicator if it is track or not
    3 recHitTheta,
    4 recHitR, 
    5 recHitX,
    6 recHitY,
    7 recHitZ,
    8 recHitTime
    (https://github.com/cms-pepr/HGCalML/blob/master/modules/datastructures/TrainData_NanoML.py#L211-L221)

    Args:
        flip (bool): If True, flips the negative endcap z-values to positive
        reduce_noise (float): Randomly delete a fraction of noise. Useful
            to speed up training.
    """
    def __init__(self, path, flip=True, reduce_noise: float=None):
        super(SinglePhotonDataset, self).__init__(path)
        self.npzs = list(sorted(glob.iglob(path + '/*.npz')))
        self.flip = flip
        self.reduce_noise = reduce_noise
        self.noise_index = -1
        self.noise_mask_cache = {}

    def blacklist(self, npzs):
        """
        Remove a list of npzs from the dataset
        Useful to remove bad events
        """
        for npz in npzs: self.npzs.remove(npz)

    def get(self, i):
        d = np.load(self.npzs[i])
        x = d['recHitFeatures']
        y = d['recHitTruthClusterIdx'].squeeze()
        if self.flip and np.mean(x[:,7]) < 0:
            # Negative endcap: Flip z-dependent coordinates
            x[:,1] *= -1 # eta
            x[:,7] *= -1 # z
        if self.reduce_noise:
            # Throw away a fraction of noise
            # Have to be careful to throw away to same noise upon
            # future calls of this function.
            mask = self.noise_mask_cache.setdefault(i, mask_fraction_of_noise(y, self.reduce_noise, self.noise_index))
            x = x[mask]
            y = y[mask]
        cluster_index = incremental_cluster_index_np(y.squeeze(), noise_index=self.noise_index)
        if np.all(cluster_index == 0): print('WARNING: No objects in', self.npzs[i])
        truth_cluster_props = np.hstack((
            d['recHitTruthEnergy'],
            d['recHitTruthPosition'],
            d['recHitTruthTime'],
            d['recHitTruthID'],
            ))
        if self.reduce_noise: truth_cluster_props = truth_cluster_props[mask]
        assert truth_cluster_props.shape == (x.shape[0], 5)
        order = cluster_index.argsort()
        return Data(
            x = torch.from_numpy(x[order]).type(torch.float),
            y = torch.from_numpy(cluster_index[order]).type(torch.int),
            truth_cluster_props = torch.from_numpy(truth_cluster_props[order]).type(torch.float),
            inpz = torch.Tensor([i])
            )

    def __len__(self):
        return len(self.npzs)
    def len(self):
        return len(self.npzs)

    def split(self, fraction):
        """
        Creates two new instances of SinglePhotonDataset with a fraction of events split
        """
        left = self.__class__(self.root, self.flip, self.reduce_noise)
        right = self.__class__(self.root, self.flip, self.reduce_noise)
        split_index = int(fraction*len(self))
        left.npzs = self.npzs[:split_index]
        right.npzs = self.npzs[split_index:]
        return left, right
