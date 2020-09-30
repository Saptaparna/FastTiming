import uproot
import awkward
import numpy as np
from scipy.sparse import csr_matrix, find
from scipy.spatial import cKDTree
from tqdm import tqdm_notebook as tqdm


import os
import os.path as osp

#print(os.environ['GNN_TRAINING_DATA_ROOT'])

#fname = 'hgcalNtuple_Sep29_RecHitEnCut.root'
fname = 'hgcalNtuple_Aug26.root'
#fname = '/home/sameasy2006/DATA/HGCNTUP_photon_3to1000_50k.root'

print(type(fname))

test = uproot.open(fname)['ana']['hgc']

#%load_ext autoreload
#%autoreload 2
import torch

from scipy.sparse import coo_matrix # to encode the cluster mappings
from sklearn.neighbors import NearestNeighbors
from datasets.graph import Graph
from datasets.graph import graph_to_sparse, save_graph

sim_indices = awkward.fromiter(test['simcluster_hits_indices'].array())
sim_indices = sim_indices[sim_indices > -1].compact()
print ('size of sim_indices = ', len(sim_indices))

sim_energy = test['simcluster_energy'].array()
sim_eta = test['simcluster_eta'].array()
sim_phi = test['simcluster_phi'].array()
#sim_Energy = test['simcluster_simEnergy'].array()
sim_pid = test['simcluster_pid'].array()
sim_xcoordinate = test['simcluster_x_impactPoint'].array()
sim_ycoordinate = test['simcluster_y_impactPoint'].array()
sim_zcoordinate = test['simcluster_z_impactPoint'].array()
sim_timecoordinate = test['simcluster_time_impactPoint'].array()

print ('simcluster_x_impactPoint = ', test['simcluster_x_impactPoint'].array())
print ('simcluster_y_impactPoint = ', test['simcluster_y_impactPoint'].array())
print ('simcluster_z_impactPoint = ', test['simcluster_z_impactPoint'].array())
print ('simcluster_time_impactPoint = ', test['simcluster_time_impactPoint'].array())

rechit_layer = test['rechit_layer'].array()
rechit_time = test['rechit_time'].array()
rechit_energy = test['rechit_energy'].array()

print ('size of rechit_energy = ', len(rechit_energy))
#print ('rechit size per event = ', rechit_energy[0].shape)
rechit_x = test['rechit_x'].array()
rechit_y = test['rechit_y'].array()
rechit_z = test['rechit_z'].array()
#print('rechit_energy = ', rechit_energy)
for i in range(0, 5):
    print ('rechit size per event = ', rechit_energy[i].shape)
    #v = (rechit_energy[i] > 0.05)
    #if(v.all()): print('rechit_energy[i] = ', rechit_energy[i])            
#rechit['rechit_layer'].content[rechit['rechit_layer'].content < 0] *= -9
rechit_x.content[rechit_z.content < 0] *= -1

evt = 1
print(len(sim_indices[evt]))
print(sim_pid[evt])
print(sim_energy[evt])
print(sim_eta[evt])
print(sim_phi[evt])
print(sim_pid[evt].size)

outbase = fname.split('/')[-1].replace('.root','')
outdir = "/".join(fname.split('/')[:-2]) + outbase + "/raw"
#outdir = "/".join(fname.split('/')[:-2]) + "/npz_hgcal_k8/" + outbase+"/raw"

skipevt = 0
for i in tqdm(range(rechit_z.size),desc='events processed'): #

    okpos = False
    okneg = False
    if (sim_eta[i][sim_eta[i] > 0]).size <=2 :
        if (sim_pid[i][sim_eta[i] > 0]).size ==2 :
            if np.unique(abs(sim_pid[i][sim_eta[i] > 0])).size == 1:
                okpos = True
        if (sim_pid[i][sim_eta[i] > 0]).size ==1 :
            if (sim_pid[i][sim_eta[i] > 0])[0] == 22:
                okpos = True

    if (sim_eta[i][sim_eta[i] < 0]).size <=2 :
        if (sim_pid[i][sim_eta[i] < 0]).size ==2 :
            if np.unique(abs(sim_pid[i][sim_eta[i] < 0])).size == 1:
                okneg = True
        if (sim_pid[i][sim_eta[i] < 0]).size ==1 :
            if (sim_pid[i][sim_eta[i] < 0])[0] == 22:
                okneg = True
    #print ('rechit_energy[i].all() before = ', rechit_energy[i].all()) 
    #print ('i = ', i)
    #print ('rechit_energy[i][i] before = ', rechit_energy[i][i])            
    #if(rechit_energy[0][0] > 0.05): print ('rechit_energy[i][0] before = ', rechit_energy[i][0]) 
    #print (type(rechit_energy))
    if (okpos | okneg) :
        #print(okpos,okneg)
        skipevt += 1

        print(sim_energy[i])
        print(sim_eta[i])
        #print ('rechit_energy[i] = ', rechit_energy[i])
#        siminfo =  np.stack((sim_energy[i],sim_eta[i],sim_phi[i]))
#        pos_siminfo = np.sum(siminfo[:,np.argwhere(sim_eta[i]>0)].squeeze(), axis=1)
#        neg_siminfo = np.sum(siminfo[:,np.argwhere(sim_eta[i]<0)].squeeze(),axis=1)
         #pos_siminfo = np.stack((np.sum(sim_energy[i][sim_eta[i]>0]),np.sum(sim_eta[i][sim_eta[i]>0])/2,np.sum(sim_phi[i][sim_eta[i]>0])/2))
         #pos_siminfo = np.stack((np.sum(sim_energy[i][sim_eta[i]>0]),np.sum(sim_eta[i][sim_eta[i]>0])/2,np.sum(sim_phi[i][sim_eta[i]>0])/2))
        #print ('sim_xcoordinate[i] = ', sim_xcoordinate[i]) 
        #print ('np.stack((np.sum(sim_xcoordinate[i])) = ', np.stack((np.sum(sim_xcoordinate[i]),np.sum(sim_ycoordinate[i][sim_eta[i]>0]))))  
        pos_siminfo = np.stack((np.sum(sim_xcoordinate[i][sim_eta[i]>0])/2,np.sum(sim_ycoordinate[i][sim_eta[i]>0])/2,np.sum(sim_zcoordinate[i][sim_eta[i]>0]/2)))
        neg_siminfo = np.stack((np.sum(sim_xcoordinate[i][sim_eta[i]<0])/2,np.sum(sim_ycoordinate[i][sim_eta[i]<0])/2,np.sum(sim_zcoordinate[i][sim_eta[i]<0]/2)))
        #neg_siminfo = np.sum(siminfo[:,np.argwhere(sim_eta[i]<0)].squeeze(),axis=1)

        print ('pos_siminfo line 106 = ', pos_siminfo)
        #print (pos_siminfo)
        break
#for i in tqdm(range(1000),desc='events processed'): #

#    if sim_pid[i].size <= 4 :
#        print('event====')
#        print(sim_pid[i])
#        print(sim_eta[i])
#        print(sim_eta[i]*sim_pid[i])
#        print(sim_energy[i])
#        print('posen:',np.sum(sim_energy[i][ sim_eta[i] > 0]))
#        print ('negen:',np.sum(sim_energy[i][ sim_eta[i] < 0]))

#        count += 1
#        if (sim_eta[i][sim_eta[i] > 0 ]).size == 2 :
#            print(sim_pid[i][sim_eta[i] > 0 ])


print (skipevt)

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

    #print(data_sel.shape)
    #print(cluster_truth.shape)


    #map to absolute index
    #print('relative indices',pairs_sel)
    pairs_sel_abs = map_idx[pairs_sel]
    #print('absolute indices',pairs_sel_abs)

    #get the types of the clusters for these edges
    incoming = cluster_truth[pairs_sel_abs[:,1],:]
    outgoing = cluster_truth[pairs_sel_abs[:,0],:]

    #print('truth shape',incoming.shape)
    #print('truth shape',outgoing.shape)

    #determine determine all edges where each edge
    #has the same non-zero category
    hads = (incoming == 1).multiply(outgoing == 1)
    ems = (incoming == 2).multiply(outgoing == 2)
    mips = (incoming == 3).multiply(outgoing == 3)

    #print('hads',hads.todense())
    #print('ems',ems.todense())
    #print('mips',mips.todense())

    tot = (hads + ems + mips).nonzero()

    #print('tot',np.unique(tot[1],return_counts=True))

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
for i in tqdm(range(rechit_z.size),desc='events processed'): #

#for i in tqdm(range(1000),desc='events processed'): #



    okpos = False
    okneg = False
    if (sim_eta[i][sim_eta[i] > 0]).size <=2 :
        if (sim_pid[i][sim_eta[i] > 0]).size ==2 :
            if np.unique(abs(sim_pid[i][sim_eta[i] > 0])).size == 1:
                okpos = True
        if (sim_pid[i][sim_eta[i] > 0]).size ==1 :
            if (sim_pid[i][sim_eta[i] > 0])[0] == 22:
                okpos = True

    if (sim_eta[i][sim_eta[i] < 0]).size <=2 :
        if (sim_pid[i][sim_eta[i] < 0]).size ==2 :
            if np.unique(abs(sim_pid[i][sim_eta[i] < 0])).size == 1:
                okneg = True
        if (sim_pid[i][sim_eta[i] < 0]).size ==1 :
            if (sim_pid[i][sim_eta[i] < 0])[0] == 22:
                okneg = True

    if not(okpos | okneg) :
        #print(okpos,okneg)
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
    #print(hit_to_clus)
    #print(np.unique(hit_to_clus,return_counts=True))
    # 0 = invalid edge, 1 = hadronic edge, 2 = EM edge, 3 = MIP edge
    cats_per_hit = cluster_cats[hit_to_clus]

    #print(cats_per_hit)

    #print(hits_in_clus.shape, hit_to_clus.shape, cats_per_hit.shape)

    hit_truth = np.stack((hits_in_clus, hit_to_clus, cats_per_hit)).T
    #hit_truth = hit_truth[np.argsort(hit_truth[:,0])]

    #print('raw hit truth',hit_truth)

    hits_to_clusters = csr_matrix((hit_truth[:,2], (hit_truth[:,0], hit_truth[:,1])),
                                  (rechit_z[i].size, np.max(hit_to_clus)+1))

    #print('sparse hit truth',hits_to_clusters.todense())
    v = rechit_energy[i] > 0.05
    pos_mask = (rechit_z[i] > 0) & (rechit_energy[i] > 0.05)
    neg_mask = ~(rechit_z[i] > 0) #and (v.all())
    #print('pos_mask line 275 = ', pos_mask)
    #neg_mask = ~pos_mask 

    rechit_indices = np.arange(rechit_z[i].size)

    pos_feats = get_features(i, pos_mask)
    neg_feats = get_features(i, neg_mask)

    #print(rechit_indices.shape, pos_mask.shape, neg_mask.shape)

    #print(rechit_indices, rechit_indices.shape)

    pos_indices = rechit_indices[pos_mask]
    neg_indices = rechit_indices[neg_mask]

    #print(pos_indices, pos_indices.shape)
    #print(neg_indices, neg_indices.shape)

    pos_coords = pos_feats[:,0:3]
    neg_coords = neg_feats[:,0:3]

    #siminfo =  np.stack((sim_energy[i],sim_eta[i],sim_phi[i]))
    #pos_siminfo = siminfo[:,np.argwhere(sim_eta[i]>0)].squeeze()
    #neg_siminfo = siminfo[:,np.argwhere(sim_eta[i]<0)].squeeze()
    outbase = fname.split('/')[-1].replace('.root','')
    outdir = "/".join(fname.split('/')[:-2]) + outbase + "/raw"
    #outdir = "/".join(fname.split('/')[:-2]) + "/npz_hgcal_50k_k8/" + outbase+"/raw"
    if not os.path.exists(outdir):
        os.makedirs(outdir)




    if okpos:
        pos_siminfo = np.stack((np.sum(sim_xcoordinate[i][sim_eta[i]>0])/2,
                                np.sum(sim_ycoordinate[i][sim_eta[i]>0])/2,
                                np.sum(sim_zcoordinate[i][sim_eta[i]>0])/2))
        #pos_siminfo = np.stack((np.sum(sim_xcoordinate[i]),
                                #np.sum(sim_ycoordinate[i][sim_eta[i]>0]),
                                #np.sum(sim_zcoordinate[i][sim_eta[i]>0])))
        #print ('pos_siminfo line 294 = ', pos_siminfo)
        # 0 = invalid edge, 1 = hadronic edge, 2 = EM edge, 3 = MIP edge
        pos_Ri, pos_Ro, pos_y = get_neighbours(pos_coords, pos_indices, hits_to_clusters)
        pos_graph = Graph(pos_feats, pos_Ri, pos_Ro, pos_y, simmatched = pos_siminfo)
        save_graph(pos_graph, '%s/%s_hgcal_graph_pos_evt%d.npz'%(outdir,outbase,i))

    if okneg:
        neg_siminfo = np.stack((np.sum(sim_xcoordinate[i][sim_eta[i]<0])/2,
                                np.sum(sim_ycoordinate[i][sim_eta[i]<0])/2,
                                np.sum(sim_zcoordinate[i][sim_eta[i]<0])/2))
        # 0 = invalid edge, 1 = hadronic edge, 2 = EM edge, 3 = MIP edge
        neg_Ri, neg_Ro, neg_y = get_neighbours(neg_coords, neg_indices, hits_to_clusters)
        neg_graph = Graph(neg_feats, neg_Ri, neg_Ro, neg_y, simmatched = neg_siminfo)
        save_graph(neg_graph, '%s/%s_hgcal_graph_neg_evt%d.npz'%(outdir,outbase,i))

#    print(pos_y)

#    pos_graph = Graph(pos_feats, pos_Ri, pos_Ro, pos_y, simmatched = np.array([]))
#    pos_graph = Graph(pos_feats, pos_Ri, pos_Ro, pos_y, simmatched = pos_siminfo)
    #print(np.unique(pos_y,return_counts=True))
#    neg_graph = Graph(neg_feats, neg_Ri, neg_Ro, neg_y, simmatched = np.array([]))
#    neg_graph = Graph(neg_feats, neg_Ri, neg_Ro, neg_y, simmatched = neg_siminfo)
    #print(np.unique(neg_y,return_counts=True))


#    print('pos_graph\n',pos_graph)
#    print('neg_graph\n',neg_graph)



    # for UnnormalizedEdgeNet
#    save_graph(pos_graph, '%s/%s_hgcal_graph_pos_evt%d.npz'%(outdir,outbase,i))
#    save_graph(neg_graph, '%s/%s_hgcal_graph_neg_evt%d.npz'%(outdir,outbase,i))


print ('total skipped events:',skipevt,'out ot total events:',rechit_z.size)

import matplotlib.pyplot as plt

def draw_sample(X, Ri, Ro, y, out,
                cmap='bwr_r',
                skip_false_edges=True,
                alpha_labels=False,
                sim_list=None):

    #let's draw only the non-noise edges
    out_mask = out > 0
    Ri = Ri[out_mask]
    Ro = Ro[out_mask]
    good_outs = out[out_mask]

    # Select the i/o node features for each segment
    feats_o = X[Ro]
    feats_i = X[Ri]
    # Prepare the figure
    fig, (ax0,ax1) = plt.subplots(1, 2, figsize=(20,12))
    cmap = plt.get_cmap(cmap)


    #if sim_list is None:
        # Draw the hits (layer, x, y)
    #    ax0.scatter(X[:,0], X[:,2], c='k')
    #    ax1.scatter(X[:,1], X[:,2], c='k')
    #else:
    ax0.scatter(X[:,0], X[:,2], c='k')
    ax1.scatter(X[:,1], X[:,2], c='k')

    # Draw the segments
    if out is not None:
        #t = tqdm.tqdm()
        color_map = {1: dict(c='blue'),
                     2: dict(c='red'),
                     3: dict(c='orange')}
        for j in range(good_outs.shape[0]):
            seg_args = color_map[out[j]]

            ax0.plot([feats_o[j,0], feats_i[j,0]],
                     [feats_o[j,2], feats_i[j,2]], '-', **seg_args)
            ax1.plot([feats_o[j,1], feats_i[j,1]],
                     [feats_o[j,2], feats_i[j,2]], '-', **seg_args)
    else:
        t = tqdm.tqdm(range(y.shape[0]))
        for j in t:
            if y[j]:
                seg_args = dict(c='b', alpha=0.4)
            elif not skip_false_edges:
                seg_args = dict(c='black', alpha=0.4)
            else: continue

            ax0.plot([feats_o[j,0], feats_i[j,0]],
                     [feats_o[j,2], feats_i[j,2]], '-', **seg_args)
            ax1.plot([feats_o[j,1], feats_i[j,1]],
                     [feats_o[j,2], feats_i[j,2]], '-', **seg_args)

    # Adjust axes
    ax0.set_xlabel('$x$ [cm]')
    ax1.set_xlabel('$y$ [cm]')
    ax0.set_ylabel('$layer$ [arb]')
    ax1.set_ylabel('$layer$ [arb]')
    plt.tight_layout()
    return fig;
