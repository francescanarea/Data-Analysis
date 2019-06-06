import random
from sympy import *
import matplotlib.pyplot as plt

f=open("eigenvectors.txt", "r")

ev1 = []
ev2 = []
evs = f.readlines()

v1 = Matrix([[evs[0].rstrip()], [evs[1].rstrip()]])
v2 = Matrix([[evs[2].rstrip()], [evs[3].rstrip()]])

v1x1 = v1.tolist()[0][0]
v1x2 = v1.tolist()[1][0]

v2x1 = v2.tolist()[0][0]
v2x2 = v2.tolist()[1][0]

x, y = [], []

for i in range(100):
	a1 = random.uniform(-1,1)	
	
	a2 = random.uniform(-1,1)

	x0 = a1*v1 + a2*v2

	x.append(x0.tolist()[0][0]) 
	y.append(x0.tolist()[1][0])

plt.plot(x, y, 'ro') 	

x, y = [-1*v1x1], [-1*v1x2]

print("V1: [" + str(v1x1) + ", " + str(v1x2)+ "]")
print("V2: [" + str(v2x1) + ", " + str(v2x2)+"]")

x.append(v1x1)
y.append(v1x2)

plt.plot(x,y)

x,y = [v2x1*-10],[v2x2*-10]

x.append(v2x1*10)
y.append(v2x2*10)

plt.plot(x,y)


plt.show() 
