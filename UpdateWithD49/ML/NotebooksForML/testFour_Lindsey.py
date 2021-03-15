import sys
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
np.set_printoptions(threshold=sys.maxsize)

raw_dir= '/home/sapta/hgcalNtuple_Aug26/clusters/'
fnamelist = [filepath for filepath in glob.glob(raw_dir+'data_*.pt')]
data_list = []
for i in tqdm(fnamelist):
    data_list.append(torch.load(i))
   
totalev = len(data_list)
print('total samples:',totalev)
#print('data_list[0].x = ', data_list[0].x)
#print('data_list[0].y = ', data_list[0].y)
#print('data_list[0].z = ', data_list[0].z)


import torch_geometric
ntrainbatch = 2 #1 #5 #10 #was set to 50
ntestbatch = 3 #1 #10
#trainloader = torch_geometric.data.DataLoader(data_list[:5], batch_size=ntrainbatch)
#testloader = torch_geometric.data.DataLoader(data_list[5:10], batch_size=ntestbatch)
trainloader = torch_geometric.data.DataLoader(data_list[:2000], batch_size=ntrainbatch)
testloader = torch_geometric.data.DataLoader(data_list[2000:totalev], batch_size=ntestbatch)

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
train_loss = []
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
    train_loss.append(np.mean(np.array(loss)))
    print('train loss:',np.mean(np.array(loss)))


from scipy.stats import norm
import matplotlib.mlab as mlab
import scipy.stats as scs
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
#%matplotlib inline

def gaussian(x,  mean,a, sigma):
    return a * np.exp(-((x - mean)**2 / (2 * sigma**2)))

true = []
pred = []
def evaluate(epoch):
        """"Evaluate the model"""
        model.zero_grad()
        torch.cuda.empty_cache()
        model.eval()
        loss = []
        frac = []
        diff = []
        for data in tqdm(testloader):
            data = data.to(device)
            result = model(data)
            lossc = resoloss(result, data.y)
            print ('result.item() = ', result[0].item())
            print ('data.y.item() = ', data.y[0].item())
            print ('result.size() = ', len(result.size()))
            resultAll = 0
            dataAll = 0.0
            for i in range(0, len(result.size())): 
                frac.append((result[i].item() - data.y[i].item())/data.y[i].item()) #print(result[i].item())
                diff.append(result[i].item() - data.y[i].item())
                #true.append(data.y[i].item())
                #pred.append(result[i].item())
                #print ('frac = ', frac[i])
                #resultAll += result[i].item()
                #dataAll += data.y[i].item()
            #resultAll = result
            #dataAll = data
            #print ('resultAll = ', resultAll)
            #print ('dataAll = ', dataAll)
            #print ('resultAll = ', resultAll/len(result.size()))
            #print ('dataAll = ', dataAll/len(result.size()))
                #print('result[i].item() = ', result[i].item())
                #print('data.y[i].item() = ', data.y[i].item())
                #frac.append((result[i].item() - data.y[i].item())/data.y[i].item())
            #frac.append((resultAll - dataAll)/dataAll)

            #print ('frac[0] = ', frac[0])
            #print('index = ', i)
            #print('result[i].item() = ', result[i].item())
            #print('data.y[i].item() = ', data.y[i].item())
            loss.append(lossc.item())
            true.append(data.y[i].item())
            pred.append(result[i].item())
        plt.plot(true[i], pred[i], 'o', markersize=5, color='blue')
        plt.legend()
        plt.ylabel('Prediction (RECO) x-coordinate')
        plt.xlabel('Truth (GEN) x-coordinate')
        plt.savefig("Figure_Scatter2D_GenVsReco"+str(epoch)+".png")
        plt.clf()
        print('test loss:',np.mean(np.array(loss)))
        fracarr = np.array(frac)#/ntestbatch
        diffcarr = np.array(diff)
        print('fracarr = ', fracarr) 
        print('diffcarr = ', diffcarr)

        #bin_heights, bin_borders, _ = plt.hist(diffcarr, bins=500, label='histogram')
        bin_heights, bin_borders, _ = plt.hist(diffcarr, bins=100, label='histogram')
        bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
        try:
            popt, _ = curve_fit(gaussian, bin_centers, bin_heights, p0=[0., 100., 1.], bounds = ([-np.inf,0,0],[np.inf,np.inf,np.inf]))
            x_interval_for_fit = np.linspace(bin_borders[0], bin_borders[-1], 100)
            plt.plot(x_interval_for_fit, gaussian(x_interval_for_fit, *popt), label='fit')
            plt.legend()
            plt.xlabel('pred - true [cm]')
            plt.ylabel('counts')
            plt.title(r'$\mathrm{pred - true:}\ \mu=%.3f,\ \sigma=%.3f\ \mathrm{for~epoch}\ %g$' %(popt[0], popt[2], epoch))
            plt.grid(True)
            plt.savefig("Figure_Random"+str(epoch)+".png")
            plt.clf()
#            plt.plot(true[i], pred[i], 'o', markersize=5, color='blue')
#            plt.legend()
#            plt.ylabel('Prediction (RECO) x-coordinate')
#            plt.xlabel('Truth (GEN) x-coordinate')
#            plt.savefig("Figure_Scatter2D_GenVsReco"+str(epoch)+".png")
#            plt.clf()
            #plt.show()

        except RuntimeError:
            print("Error - curve_fit failed")
            plt.xlabel('pred - true')
            plt.ylabel('counts')
            plt.title('pred - true / true fit failed for epoch %g' %(epoch))
            plt.grid(True)
            plt.savefig("Figure_Random"+str(epoch)+".png")
            plt.clf()
            #plt.show()

      
        for i in range(0, len(train_loss)): plt.plot(i+1, train_loss[i], 'o', color='black')
        plt.legend()
        plt.ylabel('train loss')
        plt.xlabel('epoch')
        plt.savefig("train_loss.png")
        plt.clf()

        print ('np.mean(np.array(loss)) = ', np.mean(np.array(loss)))
        return np.mean(np.array(loss))

checkpoint_dir = '/home/sapta/hgcalNtuple_Aug26/checkpoint'
os.makedirs(checkpoint_dir, exist_ok=True)
best_loss = 99999999
test_loss = []
for epoch in range(1, 201): #210
    print ('epoch:',epoch)
    train(epoch)
    loss_epoch = evaluate(epoch)
    test_loss.append(loss_epoch)
    checkpoint_file = 'model_epoch_%03i.pth.tar' % ( epoch )
    torch.save(dict(model=model.state_dict()),
                   os.path.join(checkpoint_dir,checkpoint_file ))
    if loss_epoch < best_loss:
        best_loss = loss_epoch
        print('new best test loss:',best_loss)
        torch.save(dict(model=model.state_dict()),
                   os.path.join(checkpoint_dir,'model_checkpoint_best.pth.tar' ))


#for i in range(0, len(true)): plt.plot(true[i], pred[i], 'o', markersize=5, color='blue')
##for i in range(0, len(true)): print('true[i] = ', true[i])
##for i in range(0, len(true)): print('pred[i] = ', pred[i])
#plt.legend()
#plt.ylabel('Prediction (RECO) x-coordinate')
#plt.xlabel('Truth (GEN) x-coordinate')
#plt.savefig("Figure_Scatter2D_GenVsReco.png")        
#plt.clf()

#print('test_loss = ', test_loss)
#for i in range(0, len(test_loss)): plt.plot(i+1, test_loss[i], 'o', color='blue')
#plt.legend()
#plt.ylabel('test loss')
#plt.xlabel('epoch')
#plt.savefig("test_loss.png")
#plt.clf()
        
