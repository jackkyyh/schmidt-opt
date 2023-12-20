import torch
from tqdm import tqdm
from torch.func import hessian, grad, jacfwd
import matplotlib.pyplot as plt

class Model(torch.nn.Module):
    def __init__(self, N):
        super().__init__()

        self.theta = torch.nn.Parameter(torch.tensor([torch.pi/4] * N))

    def forward(self, omega):
        return fidelity(self.theta, omega)

def fidelity(theta, omega):
    psi = psi_f(theta)
    loss = psi @ omega
    return -loss

def eta(theta, lam):
    return - constraints(theta) @ lam

def constraints(theta):
    return torch.concatenate([-theta, theta - torch.pi/2])

def psi_f(theta):
    sin = theta.sin()
    cos = theta.cos()
    vecs = torch.unbind(torch.stack([sin, cos], dim=1))
    
    equation = ','.join(chr(97 + i) for i in range(len(sin)))
    psi = torch.einsum(equation, vecs).flatten()
    return psi



def SGD_train(model, omega, epoch=50, lr=1e-2):
    loss_list = []
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    for _ in tqdm(range(epoch)):
        loss = fidelity(model.theta, omega)
        loss.backward()
        loss_list.append(loss.detach())
        opt.step()
        opt.zero_grad()
    return torch.tensor(loss_list)


def interior_point(model, omega, mu=0.1, eps_feas=1e-2, eps = 1e-2, beta=0.9, alpha=0.5):
    N = len(model.theta)
    lam = torch.ones(2*N)
    hes_f_func = hessian(fidelity)
    grad_f_func = grad(fidelity)
    jac_const_func = jacfwd(constraints)
    hes_const_func = hessian(constraints)
            
    while True:
        theta = model.theta
        e = eta(theta, lam)
        t = mu * 2*N / e
        print(f"eta={e}")
        print(f"t={t}")

        f = fidelity(theta, omega)
        grad_f= grad_f_func(theta, omega)
        hes_f = hes_f_func(theta, omega)

        const = constraints(theta)
        jac_const = jac_const_func(theta)
        hes_const = hes_const_func(theta)


        def r(theta, lam):
            grad_f= grad_f_func(theta, omega)
            jac_const = jac_const_func(theta)
            const = constraints(theta)
            r = torch.concat([
                grad_f + jac_const.T @ lam,
                - torch.diag(lam) @ const - 1/t
            ])
            return r

        term1 = hes_f + torch.tensordot(lam, hes_const, dims=1)
        term2 = jac_const.T
        term3 = - torch.diag(lam) @ jac_const
        term4 = - torch.diag(const)
        A = torch.concat([
            torch.concat([term1, term2], dim=1), 
            torch.concat([term3, term4], dim=1)
            ])

        rr = r(theta, lam)
        delta = torch.linalg.solve(A, -rr)

        d_lam = delta[N:]
        s = 0.99 * min(1, min((- lam / d_lam)[d_lam < 0])
                        if (d_lam < 0).any() else 1)
        while (constraints(theta + s * delta[:N]) >= 0).any():
            s *= beta
            
        while (r(theta + s * delta[:N], lam + s * delta[N:]).norm() > (1-alpha * s) * rr.norm()).any():
            s *= beta

        # print(s)
        with torch.no_grad():
            theta += delta[:N] * s
            model.theta.copy_(theta)
            lam += delta[N:] * s


        rr = r(theta, lam)
        r_dual = rr[:N]
        r_cent = rr[N:]
        e = eta(theta, lam)

        if(r_dual.norm() < eps_feas and r_cent.norm() < eps_feas and e < eps):
            return theta, fidelity(theta, omega)

        f = fidelity(theta, omega)
        print(f"fidelity={f}")

        const = constraints(theta)
        print(f"const={const}")

        rr = r(theta, lam).norm()
        print(f"r={rr}")

def plot_wavefunctions(psi_sgd, psi_int, omega):
    # N = int(torch.log(len(psi_sgd)))
    size = 1
    x = list(range(len(psi_sgd)))
    plt.plot(x, psi_sgd, label="SGD", marker="o", markersize=size)
    plt.plot(x, psi_int, label="interior", marker="o", markersize=size)
    plt.plot(x, omega, label="target", marker="o", markersize=size)
    # plt.title("Wavefunctions")
    plt.legend()

def plot_loss(loss_list_sgd, loss_list_int):
    plt.plot(loss_list_sgd, label="SGD")
    plt.plot(loss_list_int, label="interior")
    plt.legend()
