import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

input = open("./positions.txt")
L = float(input.readline())
N_particle = int(input.readline())
N_iter = int(input.readline())

fig = plt.figure()
ax = Axes3D(fig)
Xs = []
Ys = []
Zs = []
for t in range(N_iter):

    X=[]
    Y=[]
    Z=[]
    for i in range(N_particle):
        l = input.readline().split(' ')
        for i_s in range(len(l)):
            l[i_s] = float(l[i_s])
        X.append(l[0])
        Y.append(l[1])
        Z.append(l[2])
    Xs.append(X)
    Ys.append(Y)
    Zs.append(Z)
current_X = Xs[0]
current_Y = Ys[0]
current_Z = Zs[0]
ax.set_xlim(0,L)
ax.set_ylim(0,L)
ax.set_zlim(0,L)
wframe=None
for t in range(1,N_iter):
    if wframe:
        ax.collections.remove(wframe)
    wframe = ax.scatter(current_X, current_Y, current_Z, c='blue')
    current_X = Xs[t]
    current_Y = Ys[t]
    current_Z = Zs[t]
    plt.pause(.01)
plt.show()

input = open("./energy.txt")
N_iter = int(input.readline())
E = []
for t in range(N_iter):
    E.append(float(input.readline()))
E_mean = float(input.readline())
E_std_dev = float(input.readline())
plt.ylim(E_mean-2*E_std_dev, E_mean+2*E_std_dev)
plt.plot(range(len(E)), E)
plt.show()
