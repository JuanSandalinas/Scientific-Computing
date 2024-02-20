import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse import csr_matrix
#from methods import jacobi


class SimulationGrid:

    def __init__(self, N, D = 1, objects = None):
        """
        Creates a simulation grid.
        Inputs:
            - N: How many intervals divide space
            - D: Parameter , default set to 1
            - object: Objects positions if there are any, default set to none
        """
        self.D = 1
        self.dx = 1/N
        self.N = N
        self.initialize()


    def initialize(self):
        """
        Initializes matrix with all all 0 concentrations except 1 on first row
        """
        A = np.zeros((self.N+1, self.N+1))
        A[0,:] = 1
        self.A = A


    def time_independent(self,method_fun,*args,**kwargs):
        """
        Executes a time_step ing method given a function
        Inputs:
            - method_func: Method to use function
            - stop: Stopping criteria
            - store_step: Every how many steps store data
        """
        self.data = [self.A]
        A = np.copy(self.A)
        for A_t in method_fun(A,*args,**kwargs):
            self.data.append(A_t)


    
    def time_dependent(self,t,dt = 0.0001, time_list = [0,0.001,0.01,0.1,1.0]):
        """
        Does time dependent stepping scheme.
        Inputs:
            - t: total time
            - dt: Time step in seconds
            - time_list: times to store
        """
        self.data = [self.A]
        term = (dt*self.D)/(self.dx**2)

        if 4*term > 1:
            raise Exception ("Not stable system")

        C = np.copy(self.A)
        n_steps = int(t/dt)

        ## Matrix operations
        diagM1 = [term*np.ones(self.N),np.ones(self.N+1),term*np.ones(self.N)]
        diagM2 = [term*np.ones(self.N),np.full(self.N+1,-4*term),term*np.ones(self.N)]

        M1 = diags(diagM1 , [-1, 0,1], format = "csr")
        M2 = diags(diagM2, [-1,0,1], format = "csr")

        ## Periodic boundary conditions
        M2[-2,0] = 1
        M2[1,-1] = 1


        for k in range(n_steps):
            c1=  M1@C
            c2 = C@M2
            C[1:-1] = (c1+c2)[1:-1]
            if k*dt in time_list:
                print(k*dt)
                self.data.append(C)


    def c_analytical(self,t, i_max=100):
        """
        Returns analytical solution of diffusion equation
        Inputs:
            - y: y axis solution
            - t: time period
        """
        lst = []
        y = np.linspace(0,1,self.N)
        for y_i in y:
            result = 0
            for i in range(i_max):
                arg1 = math.erfc((1 - y_i + (2*i)) /(2 * np.sqrt(self.D*t)))
                arg2 = math.erfc((1 + y_i + (2*i)) /(2 * np.sqrt(self.D*t)))
                result += (arg1-arg2)
            lst.append(result)
        return np.array(lst)


 
if __name__ == "__main__":
    dif = SimulationGrid(50)
    y_input = np.arange(0,1,51)
    time_show = [0,0.001, 0.01, 0.1, 1] 
    dt = 0.00001
    dif.time_dependent(t = 1, dt = dt,time_list = time_show)
    for i,d in enumerate(dif.data):
        text = f"Time dependent - t = {time_show[i]}"
        plt.plot(d[::-1,1],label=text)
        plt.legend()
    plt.show()

    