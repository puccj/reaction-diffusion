import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt(open("data.dat", "rb"))
threshold = data.mean()
thresh = np.where(data > threshold, 255, 0)

extent = [0, 10, 0, 10]  #l,r,b,t

plt.imshow(thresh, cmap='hot', extent=extent)
#plt.imshow(thresh, cmap='RdBu', extent=extent)
#plt.colorbar()

plt.show()