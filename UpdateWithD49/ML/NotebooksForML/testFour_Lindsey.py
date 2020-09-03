import numpy as np
from scipy.sparse import csr_matrix, find
from scipy.spatial import cKDTree
from tqdm import tqdm_notebook as tqdm

from datasets.graph import draw_sample
import torch
import torch_geometric
import os
import os.path as osp

import glob
raw_dir= '/home/sapta/hgcalNtuple_Aug31_E10/clusters/'
fnamelist = [filepath for filepath in glob.glob(raw_dir+'data_*.pt')]
data_list = []
for i in tqdm(fnamelist):
    data_list.append(torch.load(i))
   
totalev = len(data_list)
print('data_list[0].y = ', data_list[0].y)
print('total samples:',totalev)
#print('data_list.y = ', data_list.y)

import torch_geometric
ntrainbatch = 2 #10 #was set to 50
ntestbatch = 1
trainloader = torch_geometric.data.DataLoader(data_list[:totalev-300], batch_size=ntrainbatch)
testloader = torch_geometric.data.DataLoader(data_list[totalev-300:totalev], batch_size=ntestbatch)

import os.path as osp

import torch
import torch.nn as nn
import torch.nn.functional as F

import torch_geometric.transforms as T
from torch_geometric.data import DataLoader
from tqdm import tqdm_notebook as tqdm


from models.DynamicReductionNetwork import DynamicReductionNetwork

class Net(nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        self.drn = DynamicReductionNetwork()

    def forward(self, data):
        logits = self.drn(data)
        #print ('self.drn(data) = ', logits)
        #print ('F.softplus(logits) = ', F.softplus(logits))
        #return self.drn(data)
        #return F.softplus(logits)
        return logits

device = torch.device('cuda')#('cuda' if torch.cuda.is_available() else 'cpu')
model = Net().to(device)
optimizer = torch.optim.AdamW(model.parameters(), lr=0.001)



def resoloss(output,truth):
    batch_size = output.size()[0]
    #print ('output.size()[0] = ', output.size()[0])
    mse = F.mse_loss(output, truth, reduction='mean')
    res = torch.sum((output-truth)**2/truth)/batch_size
    return (mse)
    #return (mse + 0.001*res)
    #return (res)




model.train()
def train(epoch):
    model.train()
    loss = []
    print ('len(tqdm(trainloader)) = ' , len(tqdm(trainloader)))
    print (tqdm(trainloader))
    for data in tqdm(trainloader):
            print('data.y = ', data.y)
            data = data.to(device)
            print (data)
            optimizer.zero_grad()
            result = model(data)
            print ('result = ', result)
            lossc = resoloss(result, data.y)
            loss.append(lossc.item())
            lossc.backward()
            optimizer.step()
    
    print('train loss:',np.mean(np.array(loss)))


from scipy.stats import norm
import matplotlib.mlab as mlab
import scipy.stats as scs
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
#%matplotlib inline

def gaussian(x,  mean,a, sigma):
    return a * np.exp(-((x - mean)**2 / (2 * sigma**2)))

def evaluate(epoch):
        """"Evaluate the model"""
        model.zero_grad()
        torch.cuda.empty_cache()
        model.eval()
        loss = []
        frac = []
        for data in tqdm(testloader):
            data = data.to(device)
            result = model(data)
            lossc = resoloss(result, data.y)
            print ('result.item() = ', result.item())
            print ('data.y.item() = ', data.y.item())
            frac.append((result.item() - data.y.item())/data.y.item())
            loss.append(lossc.item())


        print('test loss:',np.mean(np.array(loss)))
        fracarr = np.array(frac)

        bin_heights, bin_borders, _ = plt.hist(fracarr, bins=100, label='histogram')
        bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2

        try:
            popt, _ = curve_fit(gaussian, bin_centers, bin_heights, p0=[0., 100., 1.],bounds = ([-np.inf,0,0],[np.inf,np.inf,np.inf]))
            x_interval_for_fit = np.linspace(bin_borders[0], bin_borders[-1], 100)
            plt.plot(x_interval_for_fit, gaussian(x_interval_for_fit, *popt), label='fit')
            plt.legend()
            plt.xlabel('pred - true / true')
            plt.ylabel('counts')
            plt.title(r'$\mathrm{pred - true / true:}\ \mu=%.3f,\ \sigma=%.3f$' %(popt[0], popt[2]))
            plt.grid(True)
            plt.show()
            plt.savefig('test.png')

        except RuntimeError:
            print("Error - curve_fit failed")
            plt.xlabel('pred - true / true')
            plt.ylabel('counts')
            plt.title('pred - true / true fit failed')
            plt.grid(True)
            plt.show()

        print ('np.mean(np.array(loss)) = ', np.mean(np.array(loss)))
        return np.mean(np.array(loss))

checkpoint_dir = '/home/sapta/hgcalNtuple_Aug31_E10/checkpoint'
os.makedirs(checkpoint_dir, exist_ok=True)
best_loss = 99999999
for epoch in range(1, 2): #10
    print ('epoch:',epoch)
    train(epoch)
    loss_epoch = evaluate(epoch)
    checkpoint_file = 'model_epoch_%03i.pth.tar' % ( epoch )
    torch.save(dict(model=model.state_dict()),
                   os.path.join(checkpoint_dir,checkpoint_file ))
    if loss_epoch < best_loss:
        best_loss = loss_epoch
        print('new best test loss:',best_loss)
        torch.save(dict(model=model.state_dict()),
                   os.path.join(checkpoint_dir,'model_checkpoint_best.pth.tar' ))
