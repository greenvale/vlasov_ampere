import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation_lib
#print(animation_lib.writers.list())
from matplotlib.animation import FuncAnimation, FFMpegWriter

n_data = np.loadtxt("n.csv", delimiter=",")
nu_data = np.loadtxt("nu.csv", delimiter=",")
E_data = np.loadtxt("E.csv", delimiter=",")

nt = min(n_data.shape[0], nu_data.shape[0], E_data.shape[0])
nx_n = n_data.shape[1]
nx_nu = nu_data.shape[1]
nx_E = E_data.shape[1]

fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

lines = []
titles = ["Density n", "Momentum nu", "Electric Field E"]
data_list = [n_data, nu_data, E_data]
y_lims = [(data.min(), data.max()) for data in data_list]

for ax, title, data, ylim in zip(axs, titles, data_list, y_lims):
    line, = ax.plot([], [], lw=2)
    lines.append(line)
    ax.set_xlim(0, data.shape[1])
    ax.set_ylim(ylim)
    ax.set_title(title)
    ax.grid(True)

def update(frame):
    for line, data in zip(lines, data_list):
        line.set_data(np.arange(data.shape[1]), data[frame])
    return lines

ani = FuncAnimation(fig, update, frames=nt, blit=True, interval=50)

#writer = FFMpegWriter(fps=20, metadata=dict(artist='Me'), bitrate=1800)
#ani.save("pic_simulation.mp4", writer=writer)

plt.show()