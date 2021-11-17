import torch
import numpy as np
import glob
import uproot
# from torch_cmspepr.dataset import incremental_cluster_index_np

def get_model():
    from torch_cmspepr.gravnet_model import GravnetModelWithNoiseFilter
    model = GravnetModelWithNoiseFilter(input_dim=9, output_dim=6, k=50, signal_threshold=.05)
    ckpt = 'ckpt_train_taus_integrated_noise_Oct20_212115_best_397.pth.tar'
    model.load_state_dict(torch.load(ckpt, map_location=torch.device('cpu'))['model'])
    return model


def interface():
    t = uproot.open('hgcalNtuple_Nov16_2021.root')['ana']['hgc']
    arrays = t.arrays()
    print (arrays)
    with torch.no_grad():
        model = get_model() #torch.load('gravnetwithnoisefilter.model')
        print ('here after loading the model')
        model.eval()
        print ('here after evaluating the model')
        for i in range(t.num_entries):
            if i == 2: break
            #e = np.array(arrays['RecHitHGC_energy'][i])
            #x = np.array(arrays['RecHitHGC_x'][i])
            #y = np.array(arrays['RecHitHGC_y'][i])
            #z = np.array(arrays['RecHitHGC_z'][i])
            #t = np.array(arrays['RecHitHGC_time'][i])
            e = np.array(arrays['rechit_energy'][i])
            x = np.array(arrays['rechit_x'][i])
            y = np.array(arrays['rechit_y'][i])
            z = np.array(arrays['rechit_z'][i])
            t = np.array(arrays['rechit_time'][i])

            nhits = e.shape[0]

            r = np.sqrt(x**2 + y**2 + z**2)

            d = np.sqrt(x**2 + y**2)
            theta = np.arctan2(d, z)
            eta = -np.log(np.tan(theta/2.))

            # Make sure phi is within -pi..pi
            phi = np.arctan2(x, y) % (2.*np.pi)
            phi[phi > np.pi] -= 2.*np.pi

            features = np.vstack((
                e,
                eta,
                np.zeros_like(e),
                theta,
                r,
                x,
                y,
                z,
                t
                )).T

            assert features.shape == (nhits, 9)
            print ('here after defining features')
            # This should be the truth clustering
            #y = np.array(arrays['RecHitHGC_BestMergedSimClusterIdx'][i])
            # y = incremental_cluster_index_np(y, noise_index=-1)
            #assert y.shape == (nhits,)


            X = torch.Tensor(features)
            batch = torch.zeros(nhits, dtype=torch.long)
            print ("X = ", X)
            scores_noise_filter, pass_noise_filter, out_gravnet = model(X, batch)
            print ("scores_noise_filter = ", scores_noise_filter)
            print(scores_noise_filter, pass_noise_filter, out_gravnet)



if __name__ == '__main__':
    interface()
