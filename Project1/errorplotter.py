# -*- coding: utf-8 -*-
from pylab import*
import scipy
import numpy as np
import os
import re

def removeFileExtension(fileName):
    return re.sub(r"\.[a-zA-Z]+$", "", fileName)

for i in range(1,6): 
    fileName = "Result_n=10^" + str(i) + ".txt" 
    data = np.loadtxt(fileName, skiprows = 1)
    x = data[:,0]
    u = data[:,1]
    v = data[:,2]
    w = data[:,3]
    error = data[:,4]

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(x,error,'r-' , label = r'$\varepsilon$')
    xlabel(r'$x$', size=20, labelpad=5)
    ylabel(r'$\varepsilon(x)$', size=20, labelpad=5)
    title(r'Relative error, $\varepsilon$ for $n=10^%i$'%i  )
    ax.grid()
    ax.legend(loc='4')
    
    fileNameNoExtention = removeFileExtension(fileName)
    fileDir = os.getcwd()
    epsDir = os.path.join(fileDir,"Figures", "eps")
    pngDir = os.path.join(fileDir,"Figures", "png")
    if not os.path.exists(epsDir):
        os.makedirs(epsDir)
    if not os.path.exists(pngDir):
        os.makedirs(pngDir)
    
    fig.savefig(pngDir + "/" + fileNameNoExtention + "_error.png", bbox_inches='tight', format="png")
    fig.savefig(epsDir + "/" + fileNameNoExtention + "_error.eps", bbox_inches='tight', format="eps")
#show()