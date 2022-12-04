import numpy as np
import torch
from torch import nn
from torch.nn import functional as F

class DNN(nn.Module):
    ''' Autoencoder for dimensional reduction'''
    def __init__(self, cdim, chid1, gdim, ghid1, hid2):
        super(DNN, self).__init__()
        #self.qc = qc
        self.cdim = cdim
        self.cc1 = nn.Linear(cdim, chid1)
        self.cc2 = nn.Linear(chid1, hid2)

        self.gdim = gdim
        self.gc1 = nn.Linear(gdim, ghid1)
        self.gc2 = nn.Linear(ghid1, hid2)

        self.sigmoid = nn.Sigmoid()

    def encode(self, cm, gm):
        ch1 = F.relu(self.cc1(cm))
        ch2 = F.relu(self.cc2(ch1))

        gh1 = F.relu(self.gc1(gm))
        gh2 = F.relu(self.gc2(gh1))
        return ch2, gh2

    def forward(self, cm, gm):
        c_en, g_en = self.encode(cm, gm)
        ren = g_en.mm(c_en.t())
        #ren1 = ren[self.qc, ]
        return ren, c_en, g_en


def loss_function(recon_m, m, ltmg_m, alpha=0.5):
    loss_fn = torch.nn.MSELoss(reduction='mean')
    m.requires_grad = True
    BCE = (1-alpha) * loss_fn(recon_m, m)
    ret = (recon_m - m) ** 2
    ret = torch.mul(ret, ltmg_m)
    ret = torch.mean(ret)
    loss = BCE + alpha * ret
    return loss


def AE_function(cell_features, gene_features, exp, ltmg_m, DNN_epochs=1000, DNN_learning_rate=0.001, reg_alpha=0.5, seed = 1217):
     if torch.cuda.device_count() > 1:
        print("Let's use", torch.cuda.device_count(), "GPUs!")

    cell_exp = torch.tensor(cell_features)
    cell_exp = cell_exp.type(torch.FloatTensor)
    cdim = cell_exp.shape[1]
    chid1 = int(cdim / 3)

    gene_exp = torch.tensor(gene_features)
    gene_exp = gene_exp.type(torch.FloatTensor)
    gdim = gene_exp.shape[1]
    ghid1 = int(gdim / 3)

    gclabel = torch.tensor(exp)
    gclabel = gclabel.type(torch.FloatTensor)

    ltmg = torch.tensor(ltmg_m)
    ltmg = ltmg.type(torch.FloatTensor)

    hid2 = min(chid1, ghid1)

    torch.manual_seed(int(seed))
    model = DNN(cdim=cdim, chid1=chid1, gdim=gdim, ghid1=ghid1, hid2=hid2)
    
    if torch.cuda.is_available():
        model = model.cuda()
        cell_exp = cell_exp.cuda()
        gene_exp = gene_exp.cuda()
        loss_function = loss_function.cuda()
    
    optimizer = torch.optim.Adam(model.parameters(), lr=DNN_learning_rate)

    min_loss_val = 10000
    for epoch in range(int(DNN_epochs)):
        gcron, cf, gf = model(cell_exp, gene_exp)

        loss = loss_function(gcron, gclabel, ltmg, alpha=reg_alpha)
       

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        #if epoch % 200 == 0:
            #print("Epoch {:02d}, Loss {:9.4f}".format(epoch, loss.item()))

        if loss < min_loss_val:
            min_loss_val = loss.item()
            best_gcron = gcron
            best_cf = cf
            best_gf = gf

    best_gcron = np.array(best_gcron.detach()) 
    best_cf = np.array(best_cf.detach())
    best_gf = np.array(best_gf.detach())
    return best_gcron, best_cf, best_gf, min_loss_val

