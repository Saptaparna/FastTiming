import numpy as np
from scipy.sparse import csr_matrix, find
from scipy.spatial import cKDTree
from tqdm import tqdm_notebook as tqdm

from datasets.graph import draw_sample
import torch

model_fname = '/home/sapta/model_categ_pho.pth.tar'
#model_fname = '/home/sapta/hgcalNtuple_Aug26_8nn_siminfo/checkpoints/model_checkpoint_EdgeNetWithCategories_264403_5b5c05404f_sapta.best.pth.tar'
#model_fname = '/home/sapta/hgcalNtuple_Aug26_8nn_siminfo/checkpoints/model_checkpoint_EdgeNetWithCategories_264403_5b5c05404f_sapta.best.pth.tar'

from models.EdgeNetWithCategories import EdgeNetWithCategories
mdl = EdgeNetWithCategories(input_dim=5, hidden_dim=64, output_dim=4, n_iters=6).to('cuda:0')
mdl.load_state_dict(torch.load(model_fname)['model'])
mdl.eval()

from unionfind import UnionFind
def cluster_points(X,out):
    finder_had = UnionFind(X.shape[0])
    finder_pho = UnionFind(X.shape[0])
    finder_mip = UnionFind(X.shape[0])

    for i in (range(index.shape[0])):
        if out[i] == 1:
            finder_had.union(index[i,0], index[i,1])
        if out[i] == 2:
            finder_pho.union(index[i,0], index[i,1])
        if out[i] == 3:
            finder_mip.union(index[i,0], index[i,1])

    had_roots = np.array([finder_had.find(i) for i in range(X.shape[0])], dtype=np.uint32)
    pho_roots = np.array([finder_pho.find(i) for i in range(X.shape[0])], dtype=np.uint32)
    mip_roots = np.array([finder_mip.find(i) for i in range(X.shape[0])], dtype=np.uint32)


    had_clusters = np.unique(had_roots, return_inverse=True, return_counts=True)
    pho_clusters = np.unique(pho_roots, return_inverse=True, return_counts=True)
    mip_clusters = np.unique(mip_roots, return_inverse=True, return_counts=True)

    hads = had_clusters[0][np.where(had_clusters[2] > 4)]
    ems = pho_clusters[0][np.where(pho_clusters[2] > 4)]
    mips = mip_clusters[0][np.where(mip_clusters[2] > 4)]

    had_clusters_sel = {i: np.where(had_roots == had)[0] for i, had in enumerate(hads)}
    em_clusters_sel = {i: np.where(pho_roots == em)[0] for i, em in enumerate(ems)}
    mip_clusters_sel = {i: np.where(mip_roots == mip)[0] for i, mip in enumerate(mips)}

    return had_clusters_sel,em_clusters_sel,mip_clusters_sel

import glob
flist = [filepath for filepath in glob.iglob(r'/home/sapta/hgcalNtuple_Aug26/processed/data_*.pt')]
print (len(flist))

import torch_geometric
import torch
datacls = []
for filename in tqdm(flist):

    data = torch.load(filename).to('cuda:0')
    X = data.x.cpu().numpy()
    pred_edges = mdl(data).detach()
    pred_edges_np = pred_edges.cpu().numpy()
    index = data.edge_index.cpu().numpy().T
    Ro = index[:,0]
    Ri = index[:,1]
    y = data.y.cpu().numpy()
    out =np.argmax(pred_edges_np,axis=-1)
    z  = data.z.cpu().numpy()
    true_y = z[0]
    print ('true', true_y)
    #print ('len = ', len(cluster_points(X,out)[1].keys())) 
    if len(cluster_points(X,out)[1].keys())  > 1 :
        continue
    #print (true_y)
    for clus in cluster_points(X,out)[1].values(): print (X[clus])
    datacls.extend([torch_geometric.data.Data(x = torch.tensor(X[clus]), y=torch.tensor([true_y])) for clus in cluster_points(X,out)[1].values()])
    #print ('datacls', datacls)
processed_dir = '/home/sapta/hgcalNtuple_Aug26/clusters/'
import os
import os.path as osp
if not os.path.exists(processed_dir):
    os.makedirs(processed_dir)

for i in tqdm(range(len(datacls))):
    torch.save(datacls[i], osp.join(processed_dir, 'data_{}.pt'.format(i)))    
