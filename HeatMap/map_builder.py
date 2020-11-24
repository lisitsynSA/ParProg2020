import matplotlib.pyplot as plt
import numpy as np

with open('heatTable.txt', 'r') as f:
    table_size = f.readline()
    data = np.loadtxt('heatTable.txt', skiprows=1)
    plt.imshow(data, cmap='hot', interpolation='nearest')
    plt.show()
