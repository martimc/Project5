import sys
import matplotlib.pyplot as plt
import numpy as np
import sys
import matplotlib.font_manager

def readfile(file):
    u = []
    lines = file.readline()
    line = lines.split()
    for word in line:
        u.append(float(word))

    return u

file = open(sys.argv[1])
line1 = file.readline()
t,n = line1.split()
t = float(t); n = int(n)

array = []
for i in range(n):
    u = readfile(file)
    array.append(u)

u_x = array[50, :] #one x array
u_y = array[:, 50] #one y array

steps = np.linspace(0, 1, n+1)