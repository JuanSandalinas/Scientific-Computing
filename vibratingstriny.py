
""" 
Vibrating string with 1D wave equation. With L = 1 and c = 1 problem of Set I, Scientific computing I
Need to add the 0 conditions
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.sparse import diags
from scipy.sparse import csr_array

PLOT = True
ANIMATION = False

STEP_PLOT = 1
SAVE_ANIMATION = True

N = 100
TAO = 0.001
TOTAL_TIME = 100


def initial_cte(N,L = 1):
    """
    Given an initial conditions equations, it returns all the points
    Inputs:
        - N: Number intervals
    """
    x = np.linspace(0,L,N+1)

    return np.sin(2*np.pi*(x)),x


def matrix_A(N, tao = TAO, L = 1, c = 1):
    """
    Given a number of points, creates the matrix A of shape N,N
    Input:
        - N: Number of intervals
        - Tao: Number of steps
    """

    h = L/N
    r = (c*tao/h)**2
    diagonals = [np.full(N-2, r),np.full(N-1, 2 - 2*r), np.full(N-2, r)]
    A = diags(diagonals , [-1, 0, 1])
    A = A.toarray()
    
    return A



def time_stepping(N,tao = TAO, steps = TOTAL_TIME/TAO, plot = True):
    """
    Does the stepping scheme over time
    """
    y_new = np.zeros(N+1)
    y,x = initial_cte(N)
    y_old = np.copy(y)
    A = matrix_A(N)
    plt.plot(x,y, label = time_step)


    for time_step in range(steps):

        A_ = A@y[1:-1]
        b = -y_old[1:-1]
        y_new[1:-1] = (A_+b)
        y_old[1:-1] = np.copy(y[1:-1])
        y[1:-1] = np.copy(y_new[1:-1])


        if (time_step%STEP_PLOT) == 0:
            plt.plot(x,y, label = f"Time = {0.001*time_step}")
       


def step(interval):
    global y_new
    global y
    global y_old
    global x

    for _ in range(STEP_PLOT):
        b = -y_old[1:-1]
        y_new[1:-1] = A@y[1:-1] + b

        y_old[1:-1] = np.copy(y[1:-1])

        y[1:-1] = np.copy(y_new[1:-1])
    
    ax.clear()
    ax.set_ylim(-1.5,1.5)
    ax.set_title(f'Time_ {np.round(interval*STEP_PLOT*TAO,2)}')
    ax.plot(x,y)

    
if __name__ == "__main__":

    if PLOT == True:
        time_stepping(N)

    elif ANIMATION == True:

        fig, ax = plt.subplots()

        y_new = np.zeros(N+1)
        y,x = initial_cte(N)
        y_old =np.copy(y)       
        A = matrix_A(N)
        ax.plot(x,y)

        anim = animation.FuncAnimation(fig,step, frames=int(10_00), interval=0.001, blit=False)

        if SAVE_ANIMATION == True:
            anim.save('sin_wave_animation.mp4', fps=30)
            plt.close()


