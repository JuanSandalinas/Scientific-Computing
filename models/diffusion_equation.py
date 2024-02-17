import numpy as np
from scipy.sparse import diags
from scipy.sparse import csr_array


class SimulationGrid:
    def __init__(self, N):
        self.dx = 1/N
        self.dt = None
        self.initialize(N)
        self.shape = self.current.shape
    
    def initialize_cte(self,N):
        self.D = 1 
        self.dt = (self.dx**2)/(4*self.D)
        A = np.zeros((N+1, N+1))
        A[0,:] = 1
        self.current = A 

    def time_step(self,function,stop):
        for _ in range(self.n_steps):
            self.c = function(self.c,stop)

    def time_dependent():
        A = self.current
        dt = self.dt
        dx = self.dx
        term = (self.dt*self.D)/(self.dx**2)
        current_time = 0
        data = [A]
        times = np.arange(0,tmax+dt, dt)
        Ax = A.copy()
        nrows, ncols = A.shape
        sampling_timer = [None]
        sampling_times = []
        
        for iteration in times: 
            A0 = data[-1] 
            for k in range(1, nrows-1):
                for j in range(1, ncols-1):
                    c1 = A0[k-1, j]
                    c2 = A0[k+1, j]
                    c3 = A0[k,j-1]
                    c4 = A0[k,j+1]
                    c0 = A0[k, j]
                    Ax[k,j] = c0 + term*(c1 + c2 + c3 + c4 - (4*c0))
            Ax[:,0] = Ax[:,1]
            Ax[:,-1] = Ax[:,-2]
            sampling_timer.pop()
            
            if len(sampling_timer) < 1:
                data.append(Ax.copy())
                sampling_times.append(iteration)
                sampling_timer = [k for k in range(frequency)]
        
        return data, np.array(sampling_times)

        def c_analytical(y, t, D=1, imax=100):
            """
            Returns analytical solution of diffusion equation
            Inputs:
                - y: y axis solution
                - t: time steps
            """

            lst = []
            for y_i in y:
                result = 0
                for i in range(100):
                    arg1 = math.erfc((1 - y_i + (2*i)) /(2 * np.sqrt(D )* t))
                    arg2 = math.erfc((1 + y_i + (2*i)) /(2 * np.sqrt(D )* t))
                    result += (arg1-arg2)
                lst.append(result)
            return np.array(lst)

    
        
"""
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
        self.c_o = np.zeros((self.N+1,self.N+1))
        self.c_o[0][1:-1] = np.ones(self.N+1)
        self.c = self.c_o

    def time_step(self,function,stop):
        for _ in range(self.n_steps):
            self.c = function(self.c,stop)
"""
if __name__ == "__main__":
    dif = Diffusion(4,1)
    dif.time_step(gauss_seidel,0.01)
    print(dif.c)