

import numpy as np
from scipy.sparse import diags
from scipy.sparse import csr_array
from scipy.sparse import linalg
import math

import matplotlib.pyplot as plt


class SimulationGrid:


    def __init__(self, N, _type = 'Square'):
        """

        Creates a simulation grid.

        Inputs:

            - N: Square size, rectangle shape or circle diameter.
        """


        self.N = N-1
        self.initialize(_type)
        self.data = [np.copy(self.v)] #For simulations


    def initialize(self,  _type):
        """
        Initializes the v vector and M matrix in function of the type of grid
        """
        if _type == 'Square':
            self.m_square()
        elif _type == 'Rectangle':
            self.m_rectangle()
        elif _type == 'Circle':
            self.m_circle()
        self.eigen = np.linalg.eig(self.M)

    def m_square(self):
        """
        Considers the grid as a square of size N, without considering boundaries
        """
        self.v = np.zeros(self.N**2)
        d = [np.ones(self.v.shape[0]- self.N),np.ones(self.v.shape[0]-1),
                np.full(self.v.shape[0],-4),np.ones(self.v.shape[0]-1),np.ones(self.v.shape[0] - self.N)]
        self.M = diags(d , [-self.N,-1,0,1,self.N]).toarray()
        
    def m_rectangle(self):
        """
        Considers the grid as a rectangle of heigh of N and length 2*N, without considering boundaries
        """
        self.v = np.zeros(self.N*self.N*2)
        d = [np.ones(self.v.shape[0] - self.N*2),np.ones(self.v.shape[0]-1),
                np.full(self.v.shape[0],-4),np.ones(self.v.shape[0]-1),np.ones(self.v.shape[0]- self.N*2)]
        self.M = diags(d , [-self.N*2,-1,0,1,self.N*2]).toarray()
        
    
    def m_circle(self):
        """
        Considers the grid as a circle of diameter N, without considering boundaries
        """
        self.v = np.zeros(np.pi*self.N/2)
        d = [np.ones(self.v.shape[0] - self.N*2),np.ones(self.v.shape[0]-1),
                np.full(self.v.shape[0],-4),np.ones(self.v.shape[0]-1),np.ones(self.v.shape[0]- self.N*2)]
        self.M = diags(d , [-self.N*2,-1,0,1,self.N*2]).toarray()
        


    
    def animation(self,save_animation = False):
        """

        Animates the stepping scheme:

        Inputs:

            -   method: If using time_dependent or time_independent

            -   save_animation: True == it will save the animation, default is False
        """

        fig, ax = plt.subplots()       

        C = np.copy(self.data)


        C = np.copy(self.A)
        

        ax.imshow(C, cmap='hot', interpolation='nearest', extent=[0, 1, 0, 1])

        ax.set_xlabel('X')  

        ax.set_ylabel('Y')  

        ax.set_title('Time: 0 s') 
        

        anim = animation.FuncAnimation(fig,self.frame, fargs= (ax,), frames=int(n_steps), interval = 0.000000001)


        if save_animation == True:

            print("Starting ")

            anim.save('time_dependent_diffusion_animation.mp4', fps=60)
            plt.close()


    def frame(self, iteration, ax):

        C = self.data[iteration]

        ax.clear()

        ax.set_title(f'Time dependent(t={np.round(iteration*0.0001*50, 7)}) s')

        ax.imshow(C, cmap='hot', interpolation='nearest', extent=[0, 1, 0, 1])





if __name__ == "__main__":
    dif = SimulationGrid(3)
    print(dif.eigen)
