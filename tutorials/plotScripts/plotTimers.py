### Imports
import os
import pandas as pd
import matplotlib.pyplot as plt

directory = "Profiling"

def find_csv_files(directory):
    csv_files = []
    for filename in os.listdir(directory):
        if filename.startswith("timers_") and filename.endswith(".csv"):
            csv_files.append(os.path.join(directory, filename))
    return csv_files

csv_files = find_csv_files(directory)


data_frames = {}

for csv_file in csv_files:
    try:
        mpi_rank = int(csv_file.split("_")[-1].split(".")[0])
        data_frames[mpi_rank] = pd.read_csv(csv_file)
    except:
        pass

plt.figure(figsize=(12, 6))

unique_headers = set()
for df in data_frames.values():
    unique_headers.update(df.columns)

unique_headers = sorted(list(unique_headers - {"elapsedCpuTime"})) + ["elapsedCpuTime"]

import matplotlib.colors as mcolors
colors = list(mcolors.TABLEAU_COLORS)


nRanks = len(data_frames)
nTimers = len(unique_headers)

for mpi_rank, df in data_frames.items():
    for c, header in enumerate(unique_headers):
        if header in df.columns:
            plt.bar(header, df[header], edgecolor=colors[c], color='None')

plt.minorticks_on()
plt.ylabel('Time (s)')
plt.show()
