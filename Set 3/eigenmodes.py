if __name__ == "__main__":

    import numpy as np

    # import math

    # import matplotlib.pyplot as plt


class SimulationGrid:


    def __init__(self, N):
        """

        Creates a simulation grid.

        Inputs:

            - N: How many intervals divide space

            - D: Parameter , default set to 1

            - object_: Matrix that specifies sink points. 0 means the state changes, 1 means it does not and it is fixde
        """


        self.N = N
        #self.dt = dt

        self.initialize()
        # self.data = [np.copy(self.C)] #For simulations


    def initialize(self):
        """
        Initializes matrix with all 0 concentrations except 1 on first row
        object_ matrix specifies which points need to be updated
        """

        self.data = np.zeros((self.N+1, self.N+1))



    def m_generator(self):

        v = np.copy(self.data)

        for i in self.N:

            for j in self.N:

                if j == i:

                    v[i, j] = -4

                elif j == i - 1 or j == i + 1:

                    v[i, j] = 1

        self.data = v




    def boundary_rectangle(self, length = 1, width = 2):

        v = self.data

        for i in self.N:

            for j in self.N:

                if j < width or j > width or i < length or i > length:

                    v[i, j] = 0

        self.data = v



    def boundary_circle(self, radius = 1):

        v = self.data

        for i in self.N:

            for j in self.N:

                if abs(self.N / 2 - i) ** 2 + abs(self.N / 2 - j) ** 2 > radius:

                    v[i, j] = 0

        self.data = v



    def eigenvalue(self, Type):

        A = self.data

        if Type == 'Rectangle':

            np.linalg.eigh(A)

        else:

            np.linalg.eig(A)

        self.data = A





if __name__ == "__main__":
    dif = SimulationGrid(25)
    dif.m_generator()
    dif.eigenvalue('Rectangle')
