# -*- coding: utf-8 -*-
"""
Created on: Thu Oct 16 18:12:13 2014
@author: cst
mail to: csteinbach@vdvms.com
Abstract: Spline Interpolation with Cubic Parametric Spline
Description: Cubic Spline interpolation based on Algorithm given in Numerical Recipies
Source: Numerical Recipies
Version: 0.2
Changed on: Thu Oct 16 18:12:13 2014
TO DO´s: - 
         - 1st derivativ
         - Integration

Parameter List Description: t = independent variable, f = function value , both in ordered numpy array
                            input of values from list or matrix via classmethod 
Public Methods Description: - cspline(t,f) creates class object
                            - set_bc(bc1,bc2) defines 1st derivativ at the starting point and at the end point,
                              if not set a natural spline with 2nd derivativ = 0 is assumed
                            - set_values_by_list(t,f) creates instance of class with values given in alist
                            - set_values by matrix(M) creates instance of class with values given as [n,2] matrix
                            - evaluate() computes the 2nd derivates
                            - interpolate(x) returns interpolated value at x
                            - get_spline_function returns interpolation function for independent use
"""
from pylab import array, asarray, size, empty, ndarray
class cpspline(object):
    def __init__(self,xyz,bc1=array([None,None]),bcn=array([None,None])):
        self.__dim = len(xyz)
        if self.__dim < 2:
            raise ValueError("Too less Dimensions in Cubic Parametric Spline Init")
        
        self.set_values(xyz)
        self.__bc1 = bc1
        self.__bcn = bcn        
        self.parameterise()
        self.evaluate()

    def set_values(self,xyz):
        length = []
        for x in xyz:   
            if not isinstance(x, ndarray):
                if not isinstance(x, list):
                    raise TypeError("independent variable must be given as numpy array or list")
                else:
                    x = asarray(x)
            length.append(size(x))
            
        if not(length[1:] == length[:-1]):
            raise IndexError("length of variable and function array is not the same")
        self.__length = length[1]  
        self.__xyz = xyz
            
    def set_bc(self,bc1=None,bcn=None):
        self.__bc1 = bc1
        self.__bcn = bcn
    
    def __call__(self,tint):
        return self.interpolate(tint)

    def parameterise(self):
        self.__t    = empty([self.__length])
        self.__t[0] = 0.0
        xsq         = 0.0
        
        for i in range(1,self.__length):
            for x in self.__xyz:
                xsq = xsq + (x[i] - x[i-1] )**2
            self.__t[i] = self.__t[i-1] + xsq**0.5
            
    def get_arclength(self):
        return self.__t

    def evaluate(self):
        f   = empty([self.__length,self.__dim])
        df2 = empty([self.__length,self.__dim])
        dum = empty([self.__length])
        
        t = self.__t
        for i in range(self.__dim):
            f[:,i] = self.__xyz[i]
        
        n = self.__length
        
        for j in range(self.__dim):
            
            if (self.__bc1[j] == None):
                df2[0,j] = 0.0
                dum[0]   = 0.0
            else:
                df2[0,j] = 0.5
                dum[0]   = (3./(t[1]-t[0]))*((f[1]-f[0])/(t[1]-t[0])-self.__bc1[j])
            
            for i in range(1,n-1):
                sig      = (t[i]-t[i-1])/(t[i+1]-t[i-1])
                df2[i,j] = (sig - 1.)/(sig * df2[i-1,j] + 2.) 
                dum[i]   = (6.*((f[i+1,j]-f[i,j])/(t[i+1]-t[i])-(f[i,j]-f[i-1,j])/(t[i]-t[i-1]))/(t[i+1]-t[i-1])-sig*dum[i-1])/(sig * df2[i-1,j] + 2.)
#            dum[i] = dum[i]/(t[i+1]-t[i-1])-sig*dum[i-1]
#            dum[i] = dum[i]/(sig * df2[i-1] + 2.)
                     
            if (self.__bcn[j] == None):
                df2[n-1,j] = 0.0
                dum[n-1]   = 0.0
            else:
                df2[n-1,j] = 0.5
                dum[n-1]   = (3./(t[n-1]-t[n-2]))*(self.__bcn - (f[n-1,j]-f[n-2,j])/(t[n-1]-t[n-2]))
            
            for k in range(n-2,-1,-1):
                df2[k,j] = df2[k,j] * df2[k+1,j] + dum[k]
            
        self.__df2 = df2
            
    def interpolate(self,tint):
        t   = self.__t
        f   = empty([self.__length,self.__dim])
        for i in range(self.__dim):
            f[:,i] = self.__xyz[i]
        df2 = self.__df2                      
        KLO = 0
        KHI = self.__length - 1
        while (KHI - KLO > 1):
            K = (KHI + KLO)/2
            if (t[K] > tint):
                KHI = K
            else:
                KLO = K
                
        h = t[KHI] - t[KLO]
        if (h == 0.): raise ZeroDivisionError("Bad h value")
        A = (t[KHI] - tint)/h
        B = (tint - t[KLO])/h
        xint = []
        for i in range(self.__dim):
            xint.append(A*f[KLO,i]+B*f[KHI,i]+((A**3-A)*df2[KLO,i]+(B**3-B)*df2[KHI,i])*h**2/6.)

        return xint
        
#    def get_spline_function(self):
#        
#        def spline_function(tint):
#            t   = self.__t
#            f   = empty([self.__length,2])
#            f[:,0] = self.__x
#            f[:,1] = self.__y
#            df2 = self.__df2                      
#            KLO = 0
#            KHI = self.__length - 1
#            while (KHI - KLO > 1):
#                K = (KHI + KLO)/2
#                if (t[K] > tint):
#                    KHI = K
#                else:
#                    KLO = K
#                
#            h = t[KHI] - t[KLO]
#            if (h == 0.): raise ZeroDivisionError("Bad h value")
#            A = (t[KHI] - tint)/h
#            B = (tint - t[KLO])/h
#            xint = A*f[KLO,0]+B*f[KHI,0]+((A**3-A)*df2[KLO,0]+(B**3-B)*df2[KHI,0])*h**2/6.
#            yint = A*f[KLO,1]+B*f[KHI,1]+((A**3-A)*df2[KLO,1]+(B**3-B)*df2[KHI,1])*h**2/6.        
#
#            return xint,yint
#            
#        return spline_function
        
    def fpp(self):
        return self.__df2 
    
    def fp(self):
        pass

    def integrate(self):
        pass    

        
        
from pylab import linspace, sin, pi, plot, cos
def test():    
        
    t = linspace(0., 2.*pi, 16)
    x = empty(16)
    y = empty(16)    
    for i in range(16):
        x[i] = cos(t[i])
        y[i] = sin(t[i])

    fspline = cpspline([x,y],array)
#    f2spline = cpspline(x,y)
#    fspline.parameterise()
#    fspline.evaluate()
    arc = fspline.get_arclength()
    print(arc)
    arcf = linspace(0., arc[15], 101) 

    tfine = linspace(0., 2.*pi, 101)
    xfine = empty(101)
    yfine = empty(101)
    xint  = empty(101)
    yint  = empty(101)
#    sxint = empty(101)
#    syint = empty(101)    
#    sfnc  = fspline.get_spline_function()

    for i in range(101):
        xfine[i] = cos(tfine[i])
        yfine[i] = sin(tfine[i])
        xint[i],yint[i]  = fspline(arcf[i])   
#        sxint[i], syint[i] = sfnc(arcf[i])   
    
    plot(xfine,yfine,xint,yint)    
#    plot(xfine,ffine,sxint,syint) 
#    print(sfnc)
# test routine        
if __name__ == "__main__":  
    test()   