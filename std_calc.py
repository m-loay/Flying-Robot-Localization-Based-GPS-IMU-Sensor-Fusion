import numpy as np
print("Old Sigmas")
filename = 'log_data\\Graph1.txt' 
data = np.loadtxt(filename, delimiter=',', dtype='Float64', skiprows=1)
g1 = data[:,1]
print("The Sigma(Quad.GPS.X) = ",np.std(g1,dtype='Float64'))

filename = 'log_data\\Graph2.txt' 
data = np.loadtxt(filename, delimiter=',', dtype='Float64', skiprows=1)
g2 = data[:,1]
print("The Sigma(Quad.IMU.AX) = ",np.std(g2,dtype='Float64'))


