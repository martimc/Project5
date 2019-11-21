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
    n = int(o.readline())
    lines = o.readline()
    line = lines.split()

    for word in line:
        u.append(float(word))

    o.close()
    return u, n

u, n = readfile(sys.argv[1])

plt.rc('text', usetex=True)
plt.rc('font', family='Computer Modern', size=15)
x = np.linspace(0,1, n+1)

plt.plot(x, u)
plt.xlabel("position (x)")
plt.ylabel("Displacement (u)")
plt.title("Displacement as function of position")
plt.grid()
plt.tight_layout()
plt.savefig('x_Disp.pdf')
plt.close()
