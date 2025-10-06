import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
import pandas as pd

files = ["ho_electron_dens", "ho_electron_avgmom", "lo_elec"]

data_list = []
titles = ["Density n", "Avg momentum nu-bar", "Electric Field E"]

for file in files:
    df = pd.read_csv(f"./output/{file}.csv")
    data = df.to_numpy()
    data_list.append(data)

nt = min([data.shape[0] for data in data_list])

fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
lines = []

for ax, title, data in zip(axs, titles, data_list):
    line, = ax.plot([], [], lw=2)
    lines.append(line)
    ax.set_xlim(0, data.shape[1])
    ax.set_ylim(data.min(), data.max())
    ax.set_title(title)
    ax.grid(True)

def update(frame):
    for line, data in zip(lines, data_list):
        line.set_data(np.arange(data.shape[1]), data[frame])
    return lines

ani = FuncAnimation(fig, update, frames=nt, blit=True, interval=50)

# Save animation as MP4
#writer = FFMpegWriter(fps=20, metadata=dict(artist='Me'), bitrate=1800)
#ani.save("electron_simulation.mp4", writer=writer)

plt.show()
