import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

indat = np.loadtxt('mike_krig2',skiprows = 3)

indat = indat.reshape(320,240)

plt.figure()
plt.imshow(indat,cmap=cm.Greys)
plt.savefig('figure_kriged.pdf')
plt.savefig('figure_kriged.png')
