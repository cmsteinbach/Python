# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 13:07:39 2014
TO DO: - recusive calculation of matrix entries
@author: cst
"""
from pylab import poly1d, empty, solve, array, linspace, matrix
import matplotlib.pyplot as plt

def exactpolyinterp(f=empty(0),fp=empty(0),fpp=empty(0),fppp=empty(0)):
    fg   = int( f.size / 2    )
    pg   = int( fp.size / 2   )
    ppg  = int( fpp.size / 2  )
    pppg = int( fppp.size / 2 )
    
    polyg = fg + pg + ppg + pppg
    b     = empty([polyg])
    M     = matrix(empty([polyg,polyg]))
    
    for i in range(0,fg,1):
        for j in range(0,polyg,1):
            M[i,j] = f[i,0]**j
            b[i]   = f[i,1]
            
    for i in range(0,pg,1):
        for j in range(0,polyg,1):
            expo = j-1
            fac  = j
            k      = fg + i
            M[k,j] = max(0,j) * fp[i,0]**max(0,j-1)
            b[k]   = fp[i,1]
#            print("1st diff",expo,fac)
            
    for i in range(0,ppg,1):
        for j in range(0,polyg,1):
            expo = j-2
            fac  = j*(j-1)
            k      = fg + pg + i
            M[k,j] = max(0,(j-1)*j) * fpp[i,0]**max(0,j-2)
            b[k]   = fpp[i,1]
#            print("2nd diff",expo,fac)
            
    for i in range(0,pppg,1):
        for j in range(0,polyg,1):
            expo = j-3
            fac  = j*(j-1)*(j-2)
            k      = fg + pg + ppg + i
            M[i,j] = max(0,(j-2)*(j-1)*j) * fppp[i,0]**max(0,j-3)
            b[i]   = fppp[i,1]
#            print("3rd diff",expo,fac)
            
    c = solve(M,b)
    coeff = empty(polyg)
    for i in range(polyg):
        coeff[polyg - 1 - i] = c[i]
    
    return poly1d(coeff)
    
def plotpolynomial(p,x):
    Y = p(x)
    plt.plot(x, Y)
    plt.show()
    
    
    
def test():    
    f     = array([[-0.3,0.0],[0.0,0.3],[0.7,0.0]])
    fp    = array([[0.0,0.0],[0.7,-0.7]])    
    fpp   = array([[-0.3,0.0],[0.7,0.0]])
    pintp = exactpolyinterp(f,fp,fpp)    
        
    x = linspace(-0.3, 0.7, 100)
    plotpolynomial(pintp,x)  
    
    f     = array([[-0.3,-0.3],[0.0,0.0],[0.7,0.0]])
    fpp   = array([[-0.3,0.0],[0.7,0.0]])    
    pintp = exactpolyinterp(f,empty(0),fpp)    
        
    x = linspace(-0.3, 0.7, 100)
    plotpolynomial(pintp,x)        
    
# test routine        
if __name__ == "__main__":  
    test()        

                