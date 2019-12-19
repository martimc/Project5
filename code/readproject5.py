import sys
import matplotlib.pyplot as plt
import numpy as np
import sys
import matplotlib.font_manager

def readfile(file):
    u_forward = []
    u_backward = []
    u_CN = []
    o = open(file)
    n = o.readline()
    t,n = n.split()
    n = int(n)
    t = float(t)
    lines = o.readline()
    line = lines.split()
    for word in line:
        u_forward.append(float(word))
    lines = o.readline()
    line = lines.split()
    for word in line:
        u_backward.append(float(word))
    lines = o.readline()
    line = lines.split()
    for word in line:
        u_CN.append(float(word))


    o.close()
    return n, t, u_forward ,u_backward ,u_CN

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
