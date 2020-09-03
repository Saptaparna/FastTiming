import os
import os.path as osp
import math

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch_geometric.transforms as T
from torch_geometric.data import DataLoader
from torch_geometric.utils import normalized_cut
from torch_geometric.nn import (NNConv, graclus, max_pool, max_pool_x,
                                global_mean_pool)

torch.backends.cudnn.benchmark = True
torch.backends.cudnn.enabled = True

from datasets.hitgraphs import HitGraphDataset

import tqdm
import argparse

directed = False
sig_weight = 1.0
bkg_weight = 1.0
train_batch_size = 1
valid_batch_size = 1
n_epochs = 20
lr = 0.01
hidden_dim = 64
n_iters = 6

from training.gnn import GNNTrainer
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print('using device %s'%device)
import logging

path ='/home/sapta/hgcalNtuple_Aug31_E10/'
print(path)
full_dataset = HitGraphDataset(path, directed=directed, categorical=True)
fulllen = len(full_dataset)
tv_frac = 0.20
tv_num = math.ceil(fulllen*tv_frac)
splits = np.cumsum([fulllen-tv_num,0,tv_num])
print(fulllen, splits)

train_dataset = torch.utils.data.Subset(full_dataset,np.arange(start=0,stop=splits[0]))
valid_dataset = torch.utils.data.Subset(full_dataset,np.arange(start=splits[1],stop=splits[2]))
train_loader = DataLoader(train_dataset, batch_size=train_batch_size, pin_memory=True)
valid_loader = DataLoader(valid_dataset, batch_size=valid_batch_size, shuffle=False)

train_samples = len(train_dataset)
valid_samples = len(valid_dataset)

d = full_dataset
num_features = d.num_features
num_classes = d[0].y.dim() if d[0].y.dim() == 1 else d[0].y.size(1)


num_classes = 4
print ('num_classes',num_classes)
the_weights = np.array([1., 1., 1., 1.]) #[0.017, 1., 1., 10.]

trainer = GNNTrainer(category_weights = the_weights,
                     output_dir='/home/sapta/hgcalNtuple_Aug31_E10_8nn_siminfo', device=device)


trainer.logger.setLevel(logging.DEBUG)
strmH = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
strmH.setFormatter(formatter)
trainer.logger.addHandler(strmH)

#example lr scheduling definition
def lr_scaling(optimizer):
    from torch.optim.lr_scheduler import ReduceLROnPlateau
    return ReduceLROnPlateau(optimizer, mode='min', verbose=True,
                             min_lr=5e-7, factor=0.2,
                             threshold=0.01, patience=5)


trainer.build_model(name='EdgeNetWithCategories', loss_func='nll_loss',
                    optimizer='AdamW', learning_rate=0.001, lr_scaling=lr_scaling,
                    input_dim=num_features, hidden_dim=64, n_iters=6,
                    output_dim=num_classes)

trainer.print_model_summary()

train_summary = trainer.train(train_loader, n_epochs, valid_data_loader=valid_loader)

print(train_summary)
