# -*- coding: utf-8 -*-
from pylab import*
#from  matplotlib import pyplot as plt
import scipy
import numpy as np
import os
import re

def removeFileExtension(fileName):
    return re.sub(r"\.[a-zA-Z]+$", "", fileName)

def sourceFunction(x):
    return 100*exp(-10*x);
def secondDerivate(x):
    return 1.0 - (1 - exp(-10)) * x - exp(-10 * x);
"""
x = []
u = []
v = []

with open("Result_n=10^1.txt") as inf:
    data = inf.readlines()[1:]
    for line in data:
        parts = line.split() # split line into parts
        if len(parts) > 1:   # if at least 2 parts/columns
            # print(parts[1])   # print column 2
            x.append(parts[0])
            u.append(parts[1])
            v.append(parts[2])
"""
#my_path = os.getcwd()

for i in range(1,3):
    fileName = "Result_n=10^" + str(i) + ".txt" 
    data = np.loadtxt(fileName, skiprows = 1)
    x = data[:,0]
    u = data[:,1]
    v = data[:,2]
    w = data[:,3]
    error = data[:,4]

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    ax.plot(x,u,'b-', label = r'Analytical solution, $u(x)$ ')
    ax.plot(x,v,'r-', marker = '1' , label = r'Numerical solution, $v(x)$ ')
    xlabel(r'$x$', size=20, labelpad=5)
    ylabel(r'$f(x)$', size=20, labelpad=5)
    title(r'Aproximation for $n=10  $')
    ax.grid()
    ax.legend(loc='smart')
    fileNameNoExtention = removeFileExtension(fileName)
    fig.savefig(fileNameNoExtention + ".png", bbox_inches='tight', format="png")
    #fig.savefig('/Images/eps/' + fileNameNoExtention + ".eps", bbox_inches='tight', format="eps")
#show()