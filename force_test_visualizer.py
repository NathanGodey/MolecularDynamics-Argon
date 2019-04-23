import matplotlib.pyplot as plt
import numpy as np

input = open("./force.txt")

x = []
y = []
res = input.readline().split(' ')

while res!=['']:
    x.append(float(res[0]))
    if "e" in res[1]:
        a,b =res[1][:-2].split('e')
        y.append(float(a)*10**(int(b)))
    else:
        y.append(float(res[1]))
    res = input.readline().split(' ')
plt.plot(x,np.sqrt(y))
plt.show()
