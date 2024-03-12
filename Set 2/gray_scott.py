
import numpy as np


import matplotlib.pyplot as plt

import matplotlib.animation as animation
from scipy.sparse import diags

from scipy.sparse import csr_matrix

class Gray_scott():
    def __init__(self,N,size, dt= 1, dx = 1, Du = 0.16, Dv = 0.08, f =  0.035, k = 0.06):
        """
        Does gray scott scheme for th reaction U +2V ->3V and V-> P. 
        With dirichlet boundary conditions and periodic boundary also 
        Inputs:
            - N: Space discretization
            - n_steps: number of time steps in space
            - dx: difference scheme in x
            - Du: Diffusion coefficient
            - Dv: More diffusion coefficient
            - f: control rate
            - u: concentration u
            - v: concentration v
            - size: Square size for V
        """
        self.N = N
        self.dt = dt 
        self.dx = dx
        self.dy = dx
        self.Du = Du 
        self.Dv = Dv
        self.f = f
        self.k = k
        self.size = size
        self.initial_cte()
        self.data_u = [np.copy(self.U)]
        self.data_v = [np.copy(self.U)]
    
    def initial_cte(self):
        """
        Set's the initial cte
        Inputs:
            - size: int square size
        *Note: If size is bigger than matrix will give error
        """
        noise_u = np.random.normal(loc=0, scale=0.5 / 3)
        noise_v = np.random.normal(loc=0, scale=0.75 / 3)
        self.U = np.full((self.N+1,self.N+1),0.5+ noise_u)
        self.V = np.zeros((self.N+1,self.N+1))
        
        center_xy = (self.N+1)// 2
        
        up_xy = center_xy + self.size//2
        down_xy = center_xy - self.size//2

        self.V[down_xy:up_xy, down_xy:up_xy] = 0.25 + noise_v

    def matrices(self):
        """
        Creates the matrices for the computation
        """

        cx = self.dt/(self.dx**2)
        cy = self.dt/(self.dy**2)
        
        dx = [np.full(self.N, cx),np.full(self.N+1,-2*cx),np.full(self.N,cx)]
        dy = [np.full(self.N, cy),np.full(self.N+1,-2*cy),np.full(self.N,cy)]
        
        Mx = diags(dx , [-1,0,1]).toarray()
        My = diags(dy, [-1,0,1]).toarray()

        Mx[-2,0] = cx
        Mx[1,-1] = cx

        self.Mx = Mx
        self.My = My

        self.Mx_u = Mx*self.Du
        self.My_u = My*self.Du
        
        self.Mx_v = Mx*self.Dv
        self.My_v = My*self.Dv


    def simulation(self,stop = 0.000001, safe_data = 1):
        """
        Does one simulation
        Inputs:
            - stop: When to stop the simulation
            - safe_data: Interval to save data, default is one
        """
        self.matrices()

        ### Begin pre-computations
        f_k = self.f + self.k
        ### End pre-computations
        self.n_count = 0
        while True:
            self.n_count += 1
            term = self.U*self.V*self.V*self.dt

            Ub =  self.Du*(self.U @ self.Mx + self.My @ self.U) -  term + self.f*(1-self.U)*self.dt + self.U
            Vb =  self.Dv*(self.V @ self.Mx + self.My @ self.V) +  term - f_k*self.V*self.dt + self.V


            max_diff = np.max(np.abs(self.U[1:-1, :] - Ub[1:-1, :]))
        
            if  max_diff <= stop:
                self.U[1:-1, :] = np.copy(Ub[1:-1, :])
                self.V[1:-1, :] = np.copy(Vb[1:-1, :])
                break
            
            self.U[1:-1, :] = np.copy(Ub[1:-1, :])
            self.V[1:-1, :] = np.copy(Vb[1:-1, :])
            
            if self.n_count%safe_data == 0:
                self.data_u += [np.copy(self.U)]
                self.data_v += [np.copy(self.V)]

    def animation(self,save_animation = False):
        """

        Animates the stepping scheme:

        Inputs:
            -   method: If using time_dependent or time_independent

            -   method_fun: If time_independent, which method to use

            -   t: Total animation time

            -   dt: Time stepping size

            -   save_animation: True == it will save the animation, default is False
        """

        fig, axs = plt.subplots(1,2)
        axs = axs.flatten()     

        U = np.copy(self.data_u[0])
        V = np.copy(self.data_v[0])

        n_steps = len(self.data_u)
        axs[0].imshow(U, cmap='plasma', extent=[0, 1, 0, 1])
        axs[1].imshow(V, cmap= 'plasma', extent=[0, 1, 0, 1])

        for ax in axs:

            ax.set_xlabel('X')  

            ax.set_ylabel('Y')  

        axs[0].set_title('U Time: 0 s') 
        axs[1].set_title('V Time: 0 s')
        

        anim = animation.FuncAnimation(fig,self.frame, fargs= (axs,), frames=int(n_steps), interval = 0.00001)

        if save_animation == True:

            print("Starting ")

            anim.save('time_dependent_diffusion_animation.mp4', fps=60)
            plt.close()
        else:
            plt.show()


    def frame(self, iteration, axs):

        U = self.data_u[iteration]
        V = self.data_v[iteration]

        axs[0].clear()
        axs[1].clear()

        axs[0].set_title(f'U {iteration} ')
        axs[1].set_title(f'V {iteration}')

        axs[0].imshow(U, cmap='plasma', extent=[0, 1, 0, 1])
        axs[1].imshow(V, cmap='plasma', extent=[0, 1, 0, 1])

if __name__ == "__main__": 
    
    dif = Gray_scott(100,20)
    dif.simulation(safe_data=1)
    dif.animation(save_animation=False)
    
