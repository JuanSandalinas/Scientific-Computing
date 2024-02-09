
""" 
Vibrating string with 1D wave equation. With L = 1 and c = 1 problem of Set I, Scientific computing I
Need to add the 0 conditions
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.sparse import diags
from scipy.sparse import csr_array

PLOT = False
ANIMATION = True

STEP_PLOT = 500
SAVE_ANIMATION = False

N = 100
TAO = 0.001
TOTAL_TIME = 100

SHIFT = 10**4

def initial_cte(N,L = 1):
    """
    Given an initial conditions equations, it returns all the points
    Inputs:
        - N: Number intervals
    """
    x = np.linspace(0,L,N+1)

    return np.sin(2*np.pi*(x)),x


def matrix_A(N, tao = TAO, L = 1):
    """
    Given a number of points, creates the matrix A of shape N,N
    Input:
        - N: Number of intervals
        - Tao: Number of steps
    """

    h = L/N
    r = (tao**2)/(h**2)
    diagonals = [np.full(N-2, r),np.full(N-1, 2 - 2*r), np.full(N-2, r)]
    A = diags(diagonals , [-1, 0, 1])
    A = A.toarray()
    A = A
    
    return A



def time_stepping(N,tao = TAO, steps = TOTAL_TIME/TAO, plot = True):
    """
    Does the stepping scheme over time
    """
    y_ta = np.zeros(N+1)
    y_to,x = initial_cte(N)
    y_tb = y_to
    A = matrix_A(N)
    time_step = 0
    plt.plot(x,y_to, label = time_step)

    while True:
        A_ = A@y_to[1:-1]*SHIFT
        b = -y_tb[1:-1]*SHIFT
        y_ta[1:-1] = (A_+b)/SHIFT
        y_tb[1:-1] = y_to[1:-1]
        y_to[1:-1] = y_ta[1:-1]
        time_step += 1

        if (time_step%STEP_PLOT) == 0:
            plt.plot(x,y_to, label = f"Time = {0.001*time_step}")
        if time_step >= steps:
            break
    plt.legend()
    plt.show()

def step(interval):
    global y_ta
    global y_to
    global y_tb
    global x

    for _ in range(STEP_PLOT):

        A_ = A@(y_to[1:-1]*SHIFT)
        b = -y_tb[1:-1]*SHIFT
        y_ta[1:-1] = (A_+b)/SHIFT

        y_tb[1:-1] = y_to[1:-1]

        y_to[1:-1] = y_ta[1:-1]
    
    ax.clear()
    ax.set_ylim(-1.5,1.5)
    ax.set_title(f'Step {interval*STEP_PLOT*TAO}')
    ax.plot(x,y_to)

    
if __name__ == "__main__":

    if PLOT == True:
        time_stepping(N)

    elif ANIMATION == True:


        fig, ax = plt.subplots()

        y_ta = np.zeros(N+1)
        y_to,x = initial_cte(N)
        y_tb = y_to
       
        A = matrix_A(N)
        ax.plot(x,y_to)

        anim = animation.FuncAnimation(fig,step, frames=int(10_000), interval=0.001, blit=False)
        if SAVE_ANIMATION == True:
             anim.save('wave_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
        plt.show()

