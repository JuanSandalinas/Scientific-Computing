import numpy as np
N = 10

A = np.zeros((N+1, N+1))

A[0,:] = 1



object_ = np.copy(A)

object_[-1,:] = 1

for i in range(4):
    position = np.random.randint(0,N)
    object_[0,position] = 2
particles_position = np.where(object_ == 2)
print(particles_position)
print( np.column_stack(particles_position))
