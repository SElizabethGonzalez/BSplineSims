'''
uses the output of the python simulation to find polynomial fits for each dimension as a function of some arclength approximation
then use the polynomial fit to find the control points in julia using the following three lines of code:
f(x) = polynomial fit
basis = BSplineBAsis(5,0:5)
bfit = approximate(f,basis)

'''

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('forcpt.dat')

data1 = data[:-1] #all but last
data2 = data[1:] #starting from 1

distance = np.sqrt(np.sum(np.square(np.subtract(data2,data1)),1))
totallength = np.sum(distance)

ratio = 5/totallength
scaleddistance = ratio*distance



init =np.asarray([0])
distwzero = np.concatenate((init,scaleddistance),axis=0)

totaldist = []
fuck = 0

for i in range(0,len(distwzero)):
    fuck += distwzero[i]
    totaldist.append(fuck)

x = []
y = []
z = []


for i in range(0,len(data)):
    x.append(data[i][0])
    y.append(data[i][1])
    z.append(data[i][2])

# plt.plot(totaldist,x)
# plt.show()

# plt.plot(totaldist,y)
# plt.show()

# plt.plot(totaldist,z)
# plt.show()


xfit = np.polyfit(totaldist, x, 7)
print(xfit)
trendpoly = np.poly1d(xfit)
#print(trendpoly)
plt.plot(totaldist,x)
plt.plot(totaldist,trendpoly(totaldist))
plt.show()


yfit = np.polyfit(totaldist, y, 6)
print(yfit)
trendpoly = np.poly1d(yfit)
#print(trendpoly)
plt.plot(totaldist,y)
plt.plot(totaldist,trendpoly(totaldist))
plt.show()

zfit = np.polyfit(totaldist, z, 8)
print(zfit)
trendpoly = np.poly1d(zfit)
#print(trendpoly)
plt.plot(totaldist,z)
plt.plot(totaldist,trendpoly(totaldist))
plt.show()