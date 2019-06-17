import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

anim=False
if anim:
    input = open("./positions.txt")
    unit = float(input.readline())
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
        plt.pause(.001)
    plt.show()

p_anim = False
if p_anim:
    input = open("./momentum.txt")
    N_iter = int(input.readline())
    N_particle = int(input.readline())


    fig, ax = plt.subplots()
    Xs = []
    for t in range(N_iter):

        X=[]
        for i in range(N_particle):
            l = float(input.readline())
            X.append(l)
        Xs.append(X)
    current_X = Xs[0]
    wframe=None
    for t in range(1,N_iter):
        if wframe:
            _ = [b.remove() for b in wframe[-1]]
        wframe = ax.hist(current_X, color= 'b')
        current_X = Xs[t]
        plt.pause(.0001)
    plt.show()
input = open("./energy.txt")
Eunit = float(input.readline())
N_iter = int(input.readline())
E = []
for t in range(N_iter):
    E.append(float(input.readline())*Eunit)
E_mean = np.mean(E)
E_std_dev = stats.tstd(E)
plt.plot(range(len(E)), E, c='r')
plt.xlabel("t")
plt.ylabel("Energie")
plt.show()

input = open("./tempcin.txt")
N_iter = int(input.readline())
E = []

for t in range(N_iter):
    a = input.readline()
    E.append(float(a))
plt.plot(range(len(E)), E, c='b')
plt.xlabel("t")
plt.ylabel("Température cinétique")
plt.show()
