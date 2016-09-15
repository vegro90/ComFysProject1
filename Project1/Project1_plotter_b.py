# -*- coding: utf-8 -*-

import numpy as np
import os
import re
from matplotlib import*
def removeFileExtension(fileName):
    return re.sub(r"\.[a-zA-Z]+$", "", fileName)

def sourceFunction(x):
    return 100*exp(-10*x);
def secondDerivate(x):
    return 1.0 - (1 - exp(-10)) * x - exp(-10 * x);

x = []
u = []
v = []

with open("n=10.txt") as inf:
    data = inf.readlines()[1:]
    for line in data:
        parts = line.split() # split line into parts
        if len(parts) > 1:   # if at least 2 parts/columns
            # print(parts[1])   # print column 2
            x.append(parts[0])
            u.append(parts[1])
            v.append(parts[2])


fig = figure()
ax = fig.add_subplot(1,1,1)

ax.plot(x,u,'b-', label = r'Analytical solution, $u(x)$ ')
ax.plot(x,v,'r-', marker = '1' , label = r'Numerical solution, $v(x)$ ')
xlabel(r' $x$', size=20, labelpad=5)
ylabel(r'$f(x)$', size=20, labelpad=5)
title(r'Aproximation for $n=100$')
ax.grid()
legend(loc='smart')
#fileNameNoExtention = removeFileExtension(fileName)
#fig.savefig(fileNameNoExtention, bbox_inches='tight', format="png")
show()