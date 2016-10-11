# -*- coding: utf-8 -*-
"""
Created on: Thu May  7 15:53:18 2015
@author: cst
mail to: csteinbach@vdvms.com
Abstract: fft for 2pi periodic signal
Description: analytical fft or 2pi periodic signal
             calculates the fourier coefficents of array xj
             by recusion. xj is given as function of period T
             at n+1 equidistand nodes
             tj = j*t/n , j=0,1,2, ...n
Source: H. SOEDING
        MECHANIK V (SCHWINGUNGSLEHRE)
        VORLESUNGSMANUSKRIPT NR. 15, 1981
        IFS HAMBURG
        and
        ORIGINAL FORTRAN CODE BY S.KRUEGER TUHH M6
Version: 0.1
Changed on: Thu May  7 15:53:18 2015
TO DO´s: -
Parameter List Description: n: I*4  input   : no of nodes for one period
                            the period starts with first node at j = 0
                            and ends with the last node at n-1.
                            k: I*4 input   : order of the fourier coefficents
                            xj:R*4 input   : array with equidistant nodes
                            ak:R*4 output  : cosine - fourier coefficents
                            bk:R*4 output  : sine - fourier coefficents
Public Methods Description: -
"""
from pylab import sin,cos,pi, empty, zeros
def fft(n,k,xj):
    buff1 = float(k)/float(n-1)*2.*pi
    buff2 = 1./float(n-1)
    u     = .5*xj[0]
    u1    = .0
    c     = cos(buff1) 
    for j in range(n-1, 0, -1):
      u2 = u1
      u1 = u
      u  = 2.*c*u1 - u2 + xj[j]
    
    if (k == 0):
        ak = (u-u2)*0.5*buff2
    else:
        ak = (u-u2)*buff2
    
#    mod = n - int(n/2)*2
#    print(mod)
    mod = n%2
#    print(mod)
    if((mod > 0)and(k*2 == n)):
        ak = ak*0.5
    
    bk = 2.*u1*sin(buff1)*buff2
    
    return ak, bk
    
from matplotlib import pyplot
def test():  
    x    = empty(201)
    y    = empty(201)
    yfft = zeros(201)
    for i in range(201):
        x[i] = 2*pi/200*float(i)
        y[i] = sin(x[i])

    for k in range (5):
        ak,bk = fft(200,k,y)
        print("ak",ak)
        print("bk",bk)     
        for i in range(200):
            yfft[i] = yfft[i] + ak*cos(x[i]*float(k)) + bk*sin(x[i]*float(k))
            
    pyplot.figure()
    pyplot.plot(x,y)
    pyplot.plot(x,yfft)
   
# test routine        
if __name__ == "__main__":  
    test()