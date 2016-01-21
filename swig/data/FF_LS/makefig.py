import os
import matplotlib.pyplot as plt
import numpy as np

x, y, z = np.loadtxt("/home/ikaur/firefront/swig/data/FF_LS/final_points_200.txt", usecols = (0,1,2), unpack = True)
for i in range(len(x)):
	
	plt.scatter(x, y, c=z, alpha=0.3)


plt.grid(True)
plt.legend()
plt.show()


