### IMPORTS
import sys
import re
import matplotlib.pyplot as plt
import math

### MAIN
# Set basic plot parameters
fig, (ax1) = plt.subplots(1, sharex=False)

# Font size
fs = 16

ax1.set_ylabel("Absorbed Power (W)", fontsize=fs)
ax1.set_xlabel("Time (s)", fontsize=fs)

for item in (ax1.get_xticklabels() + ax1.get_yticklabels()):
    item.set_fontsize(14)

# Post-process each log file given in command line argument
for q in range(len(sys.argv)-1) :

    filename = sys.argv[q+1]
    file = open(filename, "r")
    lines = file.readlines()

    # Initialize empty data arrays
    times = []
    powers = []

    # Initialize dummy variables
    time = 0
    power = 0    

    # Search log file line-by-line for desired data
    for line in lines :

        if "Time = " in line and not "ExecutionTime" in line:                    
            # Pick out numbers from line string and convert to floats
            time = float(re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", line)[0])

            # Append time array with new values
            times.append(time)

        if "absorbed power" in line :                    
            # Pick out numbers from line string and convert to floats
            power = float(re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", line)[0])

            # Append time array with new values
            powers.append(power)

    ax1.plot(times, powers, lw=2.5, label=filename)

# Set miscellaneous plot parameters
ax1.legend(prop={'size': 14})
ax1.grid(True)
ax1.set_xlim(xmin=0)

plt.subplots_adjust(left=0.15,
                    bottom=0.12,
                    right=0.95,
                    top=0.95,
                    wspace=0.2,
                    hspace=0.4)

# Show plot window
plt.show()
