import sys
import matplotlib.pyplot as plt
import numpy as np
import sys
import matplotlib.font_manager
from itertools import groupby
from collections import Counter

def readfile(file):
     # opens
    u = []



    o = open(file)
    readline()
    for line in o:
        line = line.split()

        u.append(float(line[0]))

    o.close()
    return u

u = readfile(sys.argv[1])

plt.rc('text', usetex=True)
plt.rc('font', family='Computer Modern', size=15)

plt.plot(x, u)
plt.xlabel("position (x)")
plt.ylabel("Displacement (y)")
plt.title("Displacement as function of position")
plt.grid()
plt.tight_layout()
plt.savefig('plots/x_Disp.pdf')
plt.close()
