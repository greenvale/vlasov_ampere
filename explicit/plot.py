import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation


for file in ["n.csv", "nu.csv", "E.csv"]:

    data = np.loadtxt(file, delimiter=",")

    if len(data.shape) == 1:
        print(data.shape)
        print(data)
    
        plt.figure
        plt.title(file)
        plt.plot(data)
        plt.show()
    else:
        print(data.shape)
        
        nt = data.shape[0]
        nx = data.shape[1]

        fig, ax = plt.subplots()
        line, = ax.plot([], [])
        ax.set_xlim(0, nx)
        ax.set_ylim(data.min(), data.max())
        ax.set_title(file)

        def update(frame):
            line.set_data(np.arange(nx), data[frame])
            return line,

        ani = FuncAnimation(fig, update, frames=nt, blit=True, interval=20)
        plt.show()
