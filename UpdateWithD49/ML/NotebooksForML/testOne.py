import uproot
import awkward
import numpy as np
from scipy.sparse import csr_matrix, find
from scipy.spatial import cKDTree
from tqdm import tqdm_notebook as tqdm
import torch
import os
import os.path as osp


fname = 'hgcalNtuple_10k_events_Updated.root'
test = uproot.open(fname)['ana']['hgc']

from scipy.sparse import coo_matrix # to encode the cluster mappings
from sklearn.neighbors import NearestNeighbors
from datasets.graph import Graph
from datasets.graph import graph_to_sparse, save_graph

sim_indices = awkward.fromiter(test['simcluster_hits_indices'].array())
sim_indices = sim_indices[sim_indices > -1].compact()

sim_energy = test['simcluster_energy'].array()
sim_eta = test['simcluster_eta'].array()
sim_phi = test['simcluster_phi'].array()
sim_pid = test['simcluster_pid'].array()
sim_xcoordinate = test['simcluster_x_impactPoint'].array()
sim_ycoordinate = test['simcluster_y_impactPoint'].array()
sim_zcoordinate = test['simcluster_z_impactPoint'].array()
#sim_timecoordinate = test['simcluster_time_impactPoint'].array()

rechit_layer = test['rechit_layer'].array()
rechit_time = test['rechit_time'].array()
rechit_energy = test['rechit_energy'].array()

rechit_x = test['rechit_x'].array()
rechit_y = test['rechit_y'].array()
rechit_z = test['rechit_z'].array()
rechit_x.content[rechit_z.content < 0] *= -1

#print (rechit_x)

def get_category(pid):
    cats = np.zeros_like(pid) # 1 are hadrons
    cats[(pid == 22) | (np.abs(pid) == 11) | (pid == 111)] = 1 # 2 are EM showers
    cats[np.abs(pid) == 13] = 2 #3 rechit_z.sizeare MIPs
    return (cats+1) # category zero are the noise hits

def get_features(ievt, mask):
    x = rechit_x[ievt][mask]
    y = rechit_y[ievt][mask]
    layer = rechit_layer[ievt][mask]
    time = rechit_time[ievt][mask]
    energy = rechit_energy[ievt][mask]
    return np.stack((x,y,layer,time,energy)).T.astype(np.float32)

def get_neighbours(coords, map_idx, cluster_truth):
    nbrs = NearestNeighbors(algorithm='kd_tree').fit(coords)
    nbrs_sm = nbrs.kneighbors_graph(coords, 8)
    nbrs_sm.setdiag(0) #remove self-loop edges
    nbrs_sm.eliminate_zeros()
    nbrs_sm = nbrs_sm + nbrs_sm.T
    pairs_sel = np.array(nbrs_sm.nonzero()).T
    data_sel = np.ones(pairs_sel.shape[0])

    #map to absolute index
    #print('relative indices',pairs_sel)
    pairs_sel_abs = map_idx[pairs_sel]
    #print('absolute indices',pairs_sel_abs)

    #get the types of the clusters for these edges
    incoming = cluster_truth[pairs_sel_abs[:,1],:]
    outgoing = cluster_truth[pairs_sel_abs[:,0],:]

    #determine determine all edges where each edge
    #has the same non-zero category
    hads = (incoming == 1).multiply(outgoing == 1)
    ems = (incoming == 2).multiply(outgoing == 2)
    mips = (incoming == 3).multiply(outgoing == 3)



    tot = (hads + ems + mips).nonzero()



    #prepare the input and output matrices (already need to store sparse)
    r_shape = (coords.shape[0],pairs_sel.shape[0])
    eye_edges = np.arange(pairs_sel.shape[0])

    R_i = csr_matrix((data_sel,(pairs_sel[:,1],eye_edges)), r_shape, dtype=np.uint8)
    R_o = csr_matrix((data_sel,(pairs_sel[:,0],eye_edges)), r_shape, dtype=np.uint8)

    # if you address the incoming edge by the outgoing index then the edge connects two
    # hits in the same sim-cluster
    y = np.zeros(shape=pairs_sel.shape[0], dtype=np.int8)
    truth = np.squeeze(np.asarray(incoming[tot[0],tot[1]]))
    if tot[0].size > 0 and tot[1].size > 0:
        y[tot[0]] = truth

    return R_i, R_o, y


skipevt = 0
for i in tqdm(range(rechit_z.size),desc='events processed'):

    # skip events with more than 1 sim particle per endcap used for pions
    if sim_pid[i].size > 2:
        skipevt += 1
        continue

    cluster_cats = get_category(sim_pid[i])

    sim_indices_cpt = awkward.fromiter(sim_indices[i])
    if isinstance(sim_indices_cpt, np.ndarray):
        if sim_indices_cpt.size == 0: #skip events that are all noise, they're meaningless anyway
            continue
        else:
            sim_indices_cpt = awkward.JaggedArray.fromcounts([sim_indices_cpt.size],sim_indices_cpt)
    hits_in_clus = sim_indices_cpt.flatten()
    hit_to_clus = sim_indices_cpt.parents

    # 0 = invalid edge, 1 = hadronic edge, 2 = EM edge, 3 = MIP edge
    cats_per_hit = cluster_cats[hit_to_clus]


    hit_truth = np.stack((hits_in_clus, hit_to_clus, cats_per_hit)).T

    hits_to_clusters = csr_matrix((hit_truth[:,2], (hit_truth[:,0], hit_truth[:,1])),
                                  (rechit_z[i].size, np.max(hit_to_clus)+1))

    ## separate by endcaps
    pos_mask = (rechit_z[i] > 0)
    neg_mask = ~pos_mask

    rechit_indices = np.arange(rechit_z[i].size)

    pos_feats = get_features(i, pos_mask)
    neg_feats = get_features(i, neg_mask)



    pos_indices = rechit_indices[pos_mask]
    neg_indices = rechit_indices[neg_mask]



    pos_coords = pos_feats[:,0:3]
    neg_coords = neg_feats[:,0:3]

    ## collect sim info for regressing on clusters
    siminfo =  np.stack((sim_energy[i],sim_eta[i],sim_phi[i]))
    pos_siminfo = siminfo[:,np.argwhere(sim_eta[i]>0)].squeeze()
    neg_siminfo = siminfo[:,np.argwhere(sim_eta[i]<0)].squeeze()

    # 0 = invalid edge, 1 = hadronic edge, 2 = EM edge, 3 = MIP edge
    pos_Ri, pos_Ro, pos_y = get_neighbours(pos_coords, pos_indices, hits_to_clusters)
    neg_Ri, neg_Ro, neg_y = get_neighbours(neg_coords, neg_indices, hits_to_clusters)



    #build graph namedtuples
    pos_graph = Graph(pos_feats, pos_Ri, pos_Ro, pos_y, simmatched = pos_siminfo)

    neg_graph = Graph(neg_feats, neg_Ri, neg_Ro, neg_y, simmatched = neg_siminfo)



#    print('pos_graph\n',pos_graph)
#    print('neg_graph\n',neg_graph)

    outbase = fname.split('/')[-1].replace('.root','')
    outdir = "/".join(fname.split('/')[:-2]) + outbase + "/raw"
    #outdir = "/".join(fname.split('/')[:-2]) + "/npz_hgcal_pionaddsiminfo_all_k8/" + outbase+"/raw"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # for UnnormalizedEdgeNet
    save_graph(pos_graph, '%s/%s_hgcal_graph_pos_evt%d.npz'%(outdir,outbase,i))
    save_graph(neg_graph, '%s/%s_hgcal_graph_neg_evt%d.npz'%(outdir,outbase,i))


print ('total skipped events:',skipevt,'out ot total events:',rechit_z.size)
