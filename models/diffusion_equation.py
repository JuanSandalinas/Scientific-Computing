import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.sparse import diags
from scipy.sparse import csr_matrix

if __name__ == "__main__":
    from methods import jacobi, sor, gauss_seidel


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
        C = np.copy(self.A)
        n_count = 0
        for A_t in method_fun(C,*args,**kwargs):
            self.data.append(np.copy(A_t))
    
    def time_dependent_matrix(self):
        """
        Creates the matrices for time dependent difference scheme
        """

        diagonals = [np.ones(self.N), np.ones(self.N)]
        
        M1 = diags(diagonals , [-1, 1])
        M2 = diags(diagonals, [-1,1]).toarray()

        ## Cyclic boundary conditions
        M2[-2,0] = 1
        M2[1,-1] = 1

        ## Periodic boundary conditions
        
        self.M1 = M1
        self.M2 = M2

    def time_dependent_step(self,C):
        """
        Does on step of the time dependent difference scheme
        Inputs:
            - C: Matrix C, which is A over time
        Outputs:
            - C: Matrix C after one step
        """
        c1=  self.M1@C
        c2 = C@self.M2
        C[1:-1] = (C + self.term*(c1 + c2 - 4*C))[1:-1]
        return C
    

    
    def time_dependent(self,t,dt = 0.0001, time_list = [0,0.001,0.01,0.1,1.0]):
        """
        Does time dependent stepping scheme.
        Inputs:
            - t: total time
            - dt: Time step in seconds
            - time_list: times to store
        """
        self.data = [self.A]
        self.term = (dt*self.D)/(self.dx**2)
        self.time_dependent_matrix()

        if 4*self.term > 1:
            raise Exception ("Not stable system")

        C = np.copy(self.A)
        n_steps = int(t/dt)


        for k in range(n_steps):
            C = self.time_dependent_step(C)
            if k*dt in time_list:
                self.data.append(np.copy(C))


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

    def animation(self,t = 1,dt = 0.0001,save_animation = False):
        """
        Animates the stepping scheme:
        Inputs:
            -   t: Total animation time
            -   dt: Time stepping size
            -   save_animation: True == it will save the animation, default is False
        """
        fig, ax = plt.subplots()       
        C = np.copy(self.A)
        self.time_dependent_matrix()
        self.term = (dt*self.D)/(self.dx**2)
        
        n_steps = t/dt
        C = np.copy(self.A)
        
        ax.imshow(C, cmap='hot', interpolation='nearest', extent=[0, 1, 0, 1])
        ax.set_xlabel('X')  
        ax.set_ylabel('Y')  
        ax.set_title('Time: 0 s') 
        
        anim = animation.FuncAnimation(fig,self.frame, fargs= (ax, C), frames=int(n_steps), interval = 0.00001)
        plt.show()

        if save_animation == True:
            anim.save('time_diffusion_animation.mp4', fps=30)
            plt.close()

    def frame(self, iteration, ax, C):
        C = self.time_dependent_step(C)
        ax.clear()
        ax.set_title(f'Time: {np.round(iteration*0.0001)} s')
        ax.imshow(C, cmap='hot', interpolation='nearest', extent=[0, 1, 0, 1])


 
if __name__ == "__main__":
    dif = SimulationGrid(50)
    dif.animation()

    