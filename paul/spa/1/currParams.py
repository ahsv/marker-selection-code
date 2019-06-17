import numpy as np

# small mesh
#lambMesh = list(np.linspace(0.01,0.1,num=10,endpoint=True)) 
#alphaMesh = (np.linspace(0.01, 0.1, num=10, endpoint=True))

# large mesh
#lambMesh = list(np.linspace(0.1,1.0,num=10,endpoint=True))
#alphaMesh = (np.linspace(0.05, 0.5, num=10, endpoint=True))

# example mesh
lambMesh = list(np.linspace(0.01,0.1,num=10,endpoint=True)) 
alphaMesh = [0.01]

lambs = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 8.0, 12.0, 16.0, 20.0, 24.0, 28.0, 32.0]
