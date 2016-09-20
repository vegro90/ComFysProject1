# -*- coding: utf-8 -*-
from pylab import*
import scipy
import numpy as np
import os
import re

def removeFileExtension(fileName):
    return re.sub(r"\.[a-zA-Z]+$", "", fileName)

for i in range(1,4): 
    fileName = "Result_n=10^" + str(i) + ".txt" 
    data = np.loadtxt(fileName, skiprows = 1)
    x = data[:,0]
    u = data[:,1]
    v = data[:,2]
    w = data[:,3]
    error = data[:,4]

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(x,v,'r-', marker = '1' , label = r'Numerical solution, $v(x)$ ')
    ax.plot(x,u,'b-', label = r'Analytical solution, $u(x)$ ')
    xlabel(r'$x$', size=20, labelpad=5)
    ylabel(r'$f(x)$', size=20, labelpad=5)
    title(r'Aproximation for $n=$%i'%10**i  )
    ax.grid()
    ax.legend(loc='smart')
    
    fileNameNoExtention = removeFileExtension(fileName)
    fileDir = os.getcwd()
    epsDir = os.path.join(fileDir,"Figures", "eps")
    pngDir = os.path.join(fileDir,"Figures", "png")
    if not os.path.exists(epsDir):
        os.makedirs(epsDir)
    if not os.path.exists(pngDir):
        os.makedirs(pngDir)
    
    fig.savefig(pngDir + "/" + fileNameNoExtention + ".png", bbox_inches='tight', format="png")
    fig.savefig(epsDir + "/" +fileNameNoExtention + ".eps", bbox_inches='tight', format="eps")
#show()