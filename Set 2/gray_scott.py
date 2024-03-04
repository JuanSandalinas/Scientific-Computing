
import numpy as np


import matplotlib.pyplot as plt

import matplotlib.animation as animation
from scipy.sparse import diags

from scipy.sparse import csr_matrix

class Gray_scott():
    def __init__(self,N, n_steps, dt= 1, dx = 1, Du = 0.16, Dv = 0.08, f = 0.035, k = 0.06):
        """
        Does gray scott scheme for th reaction U +2V ->3V and V-> P 
        Inputs:
            - N: Space discretization
            - n_steps: number of time steps in space
            - dx: difference scheme in x
            - Du: Diffusion coefficient
            - Dv: More diffusion coefficient
            - f: control rate
            - u: concentration u
            - v: concentration v
        """
        self.N = N
        self.n_steps = n_steps
        self.dt = dt 
        self.dx = dx 
        self.Du = Du 
        self.Dv = Dv 
        self.f = f
        self.k = k
        self.U = np.zeros((N+1,N+1))
        self.V = np.zeros((N+1,N+1))
        self.C = np.copy(self.A)
        self.C[0]  = 1
        self.C[-1] = 1
    def matrices():
        cx = self.dt*self.Du/(self.dx**2)
        cy = self.dt*self.Du/(self.dy**2)
        
        dx = np.full([np.full(self.N, cx), np.full(self.N +1,-2*cx),np.full(self.N,cx)])
        dy =  np.full([np.full(self.N, cy), np.full(self.N +1,-2*cy),np.full(self.N,cy)])
        self.Mx = diags(dx , [-1,0, 1])
        self.My = diags(dy, [-1,0, 1])


    def step(self,n_steps, safe_data = 100):
        self.matrices()
        Ub =  self.U@self.My + self.Mx@self.U -  self.U*self.V**2 + f*(1-self.U)
        Vb =  self.V@self.My + self.Mx@self.V +  self.U*self.V**2 - (f+k)*self.V
        self.U = np.copy(Ub)
        self.V = np.copy(Vb)

    