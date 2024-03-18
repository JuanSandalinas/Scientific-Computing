import numpy as np

import math

import matplotlib.pyplot as plt

import matplotlib.animation as animation

from scipy.sparse import diags

from scipy.sparse import csr_matrix



class Monte_carlo_dla():
    def __init__(self, N, n_particles, p= 1, auto = True):
        """
        Creates Monte Carlo grid for simulations.
            - N: Number of intervals
            - p: probability of sticking, default p = 1
            - self.A: Matrix with all 0, 1 where particles move, and 2 on the cluster
        """
        self.N = N
        self.A = np.zeros((self.N+3,self.N+3))
        self.p = p
        self.n_particles = n_particles
        
        if auto == True:
            position = int(N/2)
            self.dla(position)
            self.random_walk()
    
    
    def dla(self,position):
        """
        Creates a sink point for DLA.

        Inputs:

            - position: Initial position
        """

        self.A[self.N, position] = 2

    
    
    def particle_creation(self):

        """
        Creates a particle for the DL. and put is in particles
        Particles are noted as a 2
        """

        self.particle_pos = [1,np.random.randint(0,self.N)]
        self.A[self.particle_pos[0],self.particle_pos[1]] = 1
    

    def clustering(self):
        """
        Returns true or false
        """
        if self.A[self.particle_pos[0]-1, self.particle_pos[1]] == 2:
          
            return True, int(0)
        elif self.A[self.particle_pos[0]+1, self.particle_pos[1]] == 2:
            
            return True, int(1)
        elif self.A[self.particle_pos[0], self.particle_pos[1]-1] == 2:
         
            return True, int(2)
        elif self.A[self.particle_pos[0], self.particle_pos[1]+1] == 2:
        
            return True, int(3)
        else:
            return False, int(-1)

    def random_walk_step(self):
        """
        Does one iteration of a random walk.
        """
        moves = np.array([[-1,0],[1,0],[0,1],[0,-1]])
        while True:

            clustering = self.clustering()
            
            self.A[self.particle_pos[0],self.particle_pos[1]] = 0
            move = np.random.randint(0,4)
            self.particle_pos += moves[move]

            if (self.particle_pos[0]) < 1 or (self.particle_pos[0] > self.N):
                break
            
            elif (self.particle_pos[1]) < 1:
                self.particle_pos[1] = self.N
                self.A[self.particle_pos[0],self.particle_pos[1]] = 1
        
            elif (self.particle_pos[1] > self.N):
                self.particle_pos[1] = 1
                self.A[self.particle_pos[0],self.particle_pos[1]] = 1

            
            elif clustering[0] == True:
                if self.p >= np.random.uniform(0,1):
                    self.A[self.particle_pos[0],self.particle_pos[1]] = 2
                    break
                else:
                    move = np.random.choice([x for x in range(4) if x != clustering[1]])
                    self.particle_pos += moves[move]
                    self.A[self.particle_pos[0],self.particle_pos[1]] = 1
                    continue
    def random_walk(self):
        """
        Initializes the random walk simulation with n_particles
            - n_particles: Number of particles generated
            - p: Probability of clustering
        """
        for _ in range(self.n_particles):
            self.particle_creation()
            self.random_walk_step()
        print(f"Done with {self.n_particles} ")
    
    def plot(self):
        """
        Plots the current matrix A
        """
        fig, ax = plt.subplots() 

        ax.imshow(self.A[1:-1, 1:-1], cmap= plt.cm.colors.ListedColormap(['darkblue', 'white', 'yellow']), interpolation='nearest', extent=[0, 1, 0, 1])
        plt.show()


if __name__ == "__main__":
    N = 100
    n_particles_list = [20_000,40_000,60_000,80_000]
    p = 1
    monte_carlos = black_scholes = np.vectorize(Monte_carlo_dla, excluded=['N'])(N,n_particles_list)


