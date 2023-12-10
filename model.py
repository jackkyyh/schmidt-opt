import torch
from tqdm import tqdm

class Model(torch.nn.Module):
    def __init__(self, N):
        super().__init__()
        
        self.theta = torch.nn.Parameter(torch.tensor([torch.pi/4] * N))
    
    def forward(self, lam):
        phi = cal_phi(self.theta)
        loss = sum(abs(phi - lam)**2)
        # loss = -(phi @ lam)
        return loss


def cal_phi(theta):
    sin = theta.abs().sin()
    cos = theta.abs().cos()
    vecs = torch.unbind(torch.stack([sin, cos], dim=1))
    
    equation = ','.join(chr(97 + i) for i in range(len(sin)))
    phi = torch.einsum(equation, vecs).flatten()
    # phi, _ = torch.sort(phi)
    return phi


def train(model, lam, epoch=50, lr=1e-2):
    loss_list = []
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    for _ in tqdm(range(epoch)):
        loss = model(lam)
        loss.backward()
        loss_list.append(loss.detach().numpy())
        opt.step()
        opt.zero_grad()
    return loss_list