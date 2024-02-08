
""" 
Vibrating string  problem of Set I, Scientific computing I
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.sparse import diags
from scipy.sparse import csr_array

PLOT = False
ANIMATION = True
FUNCTION = 2

N = 1000
TAO = 0.001

def initial_cte(N,L = 1):
    """
    Given an initial conditions equations, it returns all the points
    Inputs:
        - N: Number of points
        - L = Length of the rope
    """
    h = L/N
    x = np.linspace(0,1,N)

    return np.sin(FUNCTION*np.pi*x), x

def matrix_A(N, tao = TAO):
    """
    Given a number of points, creates teh matrix A of shape N,N
    Input:
        - N: Number of points
        - Tao: Number of steps
    """

    h = 1/N
    diagonals = [np.ones(N-1),np.full(N, 2*h**2 -2), np.ones(N-1)]
    A = diags(diagonals , [-1, 0, 1])
    A = A.toarray()
    A = (tao**2)/(h**2)*A
    
    return A


def time_stepping_plots(N,tao = TAO, steps = 10000, plot = True):
    """
    Does the stepping scheme over time
    """
    y_tb = np.zeros(N)
    y_to,x = initial_cte(N)
    A = matrix_A(N,tao)
    time_step = 0

    while True:
        plt.plot(x,y_to)
        y_ta = A@y_to - y_tb
        y_tb = y_to
        y_to = y_ta
        time_step += 1
        if time_step >= steps:
            break

def time_step_an(interval):
    global y_tb
    global y_to
    global y_ta
    global x
    global A

    y_ta = A@y_to - y_tb
    y_tb = y_to
    y_to = y_ta
    ax.clear()
    ax.set_ylim(-1.2,1.2)
    ax.plot(x,y_to)
    
    


if __name__ == "__main__":

    if PLOT == True:
        time_stepping_plot(N)

    elif ANIMATION == True:

        fig, ax = plt.subplots()

        ## PLots
        ax.set_title(f"1-D Wave equation")

        ## intial data

        y_tb = np.zeros(N)
        y_to,x = initial_cte(N)
        y_ta = np.zeros(N)
        A = matrix_A(N)
        ax.plot(x,y_to)

        animation = animation.FuncAnimation(fig, time_step_an, frames=100, interval=200, blit=False)
        plt.show()
        
