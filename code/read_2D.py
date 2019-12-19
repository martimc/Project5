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
dx = 1/n

array = []
for i in range(n+1):
    u = readfile(file)
    array.append(u)

array = np.asarray(array)

u_x = array[50,:] #one x array
u_y = array[:,50] #one y array

steps = np.linspace(0, 1, n+1)

plt.plot(steps, u_y, label='displacement in y-direction for x = 0.5')
plt.plot(steps, u_x, label='displacement in x-direction for y = 0.5')
plt.xlabel('length')
plt.ylabel('u(x,y,t)')
plt.title('t = %.2fs, with dx = %.2f' % (t, dx))
plt.legend()
plt.savefig('2D.pdf')
plt.close()