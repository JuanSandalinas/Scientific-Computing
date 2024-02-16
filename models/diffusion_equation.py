import numpy as np
from scipy.sparse import diags
from scipy.sparse import csr_array
from methods import gauss_seidel,sor

class Diffusion():

    def __init__(self, N, T,L = 1, tao = 0.001,auto = True):
        self.N = N
        self.T = T
        self.tao = tao

        self.n_steps = int(T/tao)
        self.x = np.linspace(0,L,self.N+1)
        self.y = np.linspace(0,L, self.N+1)

        if auto == True:
            self.initial_cte()

    def initial_cte(self):
        self.c_o = np.zeros((self.N+1,self.N+3))
        self.c_o[0][1:-1] = np.ones(self.N+1)
        self.c = self.c_o

    def time_step(self,function,stop):
        for _ in range(self.n_steps):
            self.c = function(self.c,stop)

if __name__ == "__main__":
    dif = Diffusion(4,1)
    dif.time_step(gauss_seidel,0.01)
    print(dif.c)