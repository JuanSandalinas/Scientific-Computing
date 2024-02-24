
""" 
Vibrating string with 1D wave equation. With L = 1 and c = 1 problem of Set I, Scientific computing I
Need to add the 0 conditions
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.sparse import diags
from scipy.sparse import csr_array


class Vibrating_string:
    """
    Vibrating string  class:
    Inputs:
        - mode: mode == 1 uses sin(2*pi*x), mode == 2 uses sin(5*pi*x), mode == 3 uses sin
        - N: Number of intervals
        - T: Total time of the simulation in seconds (s)
    """

    def __init__(self,mode, N, T, L = 1, tao = 0.001, c = 1, auto = True):
        self.mode = mode
        
        self.N = N
        self.T = T
        self.L = L
        self.tao = tao
        self.c = c

        self.h =L/N
        self.n_steps = int(T/tao)
        self.data = []

        if auto == True:
            self.initial_cte()
            self.iteration_matrix()
            self.time_stepping()

    def initial_cte(self):
        """
        It creates the initial conditions
        """
        self.x = np.linspace(0,self.L,self.N+1)
        self.y_o = np.zeros(self.N +1)
       
        if self.mode == 1:
            self.y_o = np.sin(2*np.pi*self.x)
            self.y = self.y_o
            
        elif self.mode == 2:
            self.y_o = np.sin(5*np.pi*self.x)
            self.y = self.y_o

        elif self.mode ==3:
            pos = ((1/5) < self.x) & ((2/5)> self.x)
            self.y_o[pos] = np.sin(5 * np.pi * self.x[pos])
            self.y = self.y_o


    def iteration_matrix(self):
        """
        Creates the iteration matrix A.
        """

        r = (self.c*self.tao/self.h)**2
        diagonals = [np.full(self.N-2, r),np.full(self.N-1, 2 - 2*r), np.full(self.N-2, r)]
        self.A = diags(diagonals , [-1, 0, 1])

    def step(self):
        """ 
        Stepping scheme int two steps:
            1. Iterates in the form y(t+1) = A*y(t) - y(t-1)
            2. Moves to next time t, y(t-1) = y(t), and y(t) = y(t+1)
        Note:
            - Iteration are done without including the boundaries since they are 0
        """
        self.y_new[1:-1] = self.A@self.y[1:-1] - self.y_old[1:-1]
        self.data.append(np.copy(self.y_new))
        self.y_old[1:-1] = self.y[1:-1]
        self.y[1:-1] = self.y_new[1:-1]
    
    def time_stepping(self):
        self.y_new = np.zeros(self.N+1)
        self.y_old = np.copy(self.y_o)
        for time_step in range(self.n_steps):
            self.step()



    def animation(self,save_animation = False):
        """
        Animates the stepping scheme:
        Inputs:
            -   save_animation: True == it will save the animation, default is False
        """
        fig, ax = plt.subplots()
        self.y_new = np.zeros(self.N+1)
        self.y_old =np.copy(self.y_o)       

        ax.plot(self.x,self.y_o)

        anim = animation.FuncAnimation(fig,self.frame, fargs= (ax,), frames=int(self.n_steps), interval = 0.0000000001)
        if save_animation == True:
            anim.save('sin5pi_if_wave_animation.mp4', fps=30)
            plt.close()
        

    

    def frame(self, iteration, ax):
        ax.clear()
        ax.set_ylim(-1.5,1.5)
        ax.set_title(f'Time = {np.round(iteration*self.tao,6)}')
        ax.plot(self.x,self.data[iteration])



if __name__ == "__main__":
    string = Vibrating_string(3,100,1)
    string.animation(save_animation = True)
    


