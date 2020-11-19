import numpy as np
import matplotlib.pyplot as plt

sols = np.genfromtxt("sols.csv",delimiter=",")

fig, [ax1,ax2] = plt.subplots(nrows=1,ncols=2)

ax1.plot(sols[:,0])
ax2.plot(sols[:,0],sols[:,1])

plt.show()
