import sys
import matplotlib.pyplot as plt
import numpy as np
import sys
import matplotlib.font_manager



def f(x,t):
    N = 1000
    u = x
    for n in range(1,N):
        u += (-1)**n*2/(n*np.pi)*np.sin(n*np.pi*x)*np.exp(-1*(n*np.pi)**2*t)
    return u

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
    return n, t, np.array(u_forward) ,np.array(u_backward) ,np.array(u_CN)

n, t, u_forward, u_backward, u_CN = readfile(sys.argv[1])
x = np.linspace(0,1,n+1)

FEerror = np.max(np.abs(f(x,t) - u_forward))
IEerror = np.max(np.abs(f(x,t) - u_backward))
CNerror = np.max(np.abs(f(x,t) - u_CN))


print('The maximum error for FE: {}'.format(FEerror))
print('The maximum error for IE: {}'.format(IEerror))
print('The maximum error for CN: {}'.format(CNerror))
print('Difference between CN and FE: {}'.format(np.max(np.abs(u_CN - u_forward))))


"""
plt.rc('text', usetex=True)
plt.rc('font', family='Computer Modern', size=15)

plt.plot(x, u_forward,'b', label = 'Forward Euler')
plt.plot(x, u_backward,'r', label = 'Implicit Euler')
plt.plot(x, u_CN,'y', label = 'Crank-Nicholson')
plt.legend()
plt.xlabel("position x")
plt.ylabel("Displacement u(t)".format(t))
plt.title("Solution for 3 schemes at t = {}".format(t))
plt.grid()
plt.axis('equal')
plt.tight_layout()
plt.savefig('t50x1.pdf')
plt.close()

plt.plot(x, u_backward)
plt.xlabel("position x")
plt.ylabel("Displacement u(t)".format(t))
plt.title("Implicit Euler solution for t = {}".format(t))
plt.grid()
plt.axis('equal')
plt.tight_layout()
plt.savefig('IE.pdf')
plt.close()

plt.plot(x, u_CN)
plt.xlabel("position x")
plt.ylabel("Displacement u(t)".format(t))
plt.title("Crank-Nicholson solution for t = {}".format(t))
plt.grid()
plt.axis('equal')
plt.tight_layout()
plt.savefig('CN.pdf')
plt.close()
"""
