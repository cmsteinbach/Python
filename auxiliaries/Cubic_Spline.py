# -*- coding: utf-8 -*-
"""
Created on: Thu Oct 16 18:12:13 2014
@author: cst
mail to: csteinbach@vdvms.com
Abstract: Spline Interpolation with Cubic Spline
Description: Cubic Spline Paramerisation based on Algorithm given in Numerical Recipies
Source: Numerical Recipies
Version: 0.2
Changed on: Thu Jul 30 2015
TO DO´s: - 1st derivativ
         - Integration

Parameter List Description: t = independent variable, f = function value , both in ordered numpy array or list

Public Methods Description: - cspline(t,f) creates class object
                            - set_bc(bc1,bc2) defines 1st derivativ at the starting point and at the end point,
                              if not set a natural spline with 2nd derivativ = 0 is assumed
                            - set_values_by_list(t,f) creates instance of class with values given in alist - removed
                            - set_values by matrix(M) creates instance of class with values given as [n,2] matrix - removed
                            - evaluate() computes the 2nd derivates
                            - interpolate(x) returns interpolated value at x
                            - get_spline_function returns interpolation function for independent use - removed
"""
from pylab import asarray, size, empty, ndarray
class cspline(object):
    def __init__(self,t,f,bc1=None,bcn=None):
        self.set_values(t,f)
        self.__bc1 = bc1
        self.__bcn = bcn  
        self.evaluate()
            
    def set_bc(self,bc1=None,bcn=None):
        self.__bc1 = bc1
        self.__bcn = bcn
        
    def set_values(self,t,f):
        if not isinstance(t, ndarray):
            if not isinstance(t, list):
                raise TypeError("independent variable must be given as numpy array or list")
            else:
               self.__t = asarray(t)
        else:
            self.__t = t
        if not isinstance(f, ndarray):
            if not isinstance(f, list):
                raise TypeError("function value must be given as numpy array or list")
            else:
               self.__f = asarray(f) 
        else:
            self.__f = f
                 
        f_len = size(f)
        t_len = size(t)
        if not (f_len==t_len):
            raise IndexError("length of variable and function array is not the same")
        self.__length = t_len        
    
#    @classmethod
#    def set_values_by_list(cls,t,f):
#        if not isinstance(t, list):
#            raise TypeError("independent variable must be given as list")
#        if not isinstance(f, list):
#            raise TypeError("function value must be given as numpy list")        
#        
#        t_array = asarray(t)
#        f_array = asarray(f)
#
#        return cls(t_array,f_array)

#    @classmethod
#    def set_values_by_matrix(cls,M):
#        if not isinstance(M, ndarray):
#            raise TypeError("independent variable and function values must be given as numpy array")     
#        
#        t_array = M[:,0]
#        f_array = M[:,1]

#        return cls(t_array,f_array)

    def evaluate(self):
        df2 = empty([self.__length])
        dum = empty([self.__length])
        
        t = self.__t
        f = self.__f
        n = self.__length
        
        if (self.__bc1 == None):
            df2[0] = 0.0
            dum[0] = 0.0
        else:
            df2[0] = 0.5
            dum[0] = (3./(t[1]-t[0]))*((f[1]-f[0])/(t[1]-t[0])-self.__bc1)
            
        for i in range(1,n-1):
            sig    = (t[i]-t[i-1])/(t[i+1]-t[i-1])
            df2[i] = (sig - 1.)/(sig * df2[i-1] + 2.) 
            dum[i] = (6.*((f[i+1]-f[i])/(t[i+1]-t[i])-(f[i]-f[i-1])/(t[i]-t[i-1]))/(t[i+1]-t[i-1])-sig*dum[i-1])/(sig * df2[i-1] + 2.)
#            dum[i] = dum[i]/(t[i+1]-t[i-1])-sig*dum[i-1]
#            dum[i] = dum[i]/(sig * df2[i-1] + 2.)
                     
        if (self.__bcn == None):
            df2[n-1] = 0.0
            dum[n-1] = 0.0
        else:
            df2[n-1] = 0.5
            dum[n-1] = (3./(t[n-1]-t[n-2]))*(self.__bcn - (f[n-1]-f[n-2])/(t[n-1]-t[n-2]))
            
        for k in range(n-2,-1,-1):
            df2[k] = df2[k] * df2[k+1] + dum[k]
            
        self.__df2 = df2
            
    def interpolate(self,tint):
        KLO = 0
        KHI = self.__length - 1
        t   = self.__t
        f   = self.__f
        df2 = self.__df2
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
        fint = A*f[KLO]+B*f[KHI]+((A**3-A)*df2[KLO]+(B**3-B)*df2[KHI])*h**2/6.
        return fint

    def __call__(self,tint):
        return self.interpolate(tint)
        
#    def get_spline_function(self):
#        class spline_function(object):
#            def __init__(self,KHI,t,f,df2):
#                self.__KHI = KHI
#                self.__t   = t
#                self.__f   = f
#                self.__df2 = df2                
#            def __call__(self,tint):
#                KLO = 0
#                KHI = self.__KHI
#                t   = self.__t
#                f   = self.__f
#                df2 = self.__df2
#                while (KHI - KLO > 1):
#                    K = (KHI + KLO)/2
#                    if (t[K] > tint):
#                        KHI = K
#                    else:
#                        KLO = K
#                
#                h = t[KHI] - t[KLO]
#                if (h == 0.): raise ZeroDivisionError("Bad h value")
#                A = (t[KHI] - tint)/h
#                B = (tint - t[KLO])/h
#                fint = A*f[KLO]+B*f[KHI]+((A**3-A)*df2[KLO]+(B**3-B)*df2[KHI])*h**2/6.  
#                return fint
                
#        return spline_function(self.__length - 1,self.__t,self.__f,self.__df2)
        
    def fpp(self):
        return self.__df2 
    
    def fp(self):
        pass

    def integrate(self):
        pass    

        
        
from pylab import linspace, sin, cos, pi, plot
def test():    
        
    x = linspace(0., 2.*pi, 6)
    f = empty(6)
    for i in range(6):
        f[i] = sin(x[i])

    fspline = cspline(x,f)
    fspline.evaluate()

    xfine = linspace(0., 2.*pi, 101)
    ffine = empty(101)
    fint  = empty(101)
    sfint = empty(101)    
#    sfnc  = fspline.get_spline_function()
#    print(sfnc)

    for i in range(101):
        ffine[i] = sin(xfine[i])
        fint[i]  = fspline(xfine[i])   
#        sfint[i] = sfnc(xfine[i])   
        
#    for i in range(6):
#        f[i] = cos(x[i])
        
#    fspline.set_values(x,f)
#    fspline.evaluate()
#    for i in range(101):
#        fint[i]  = fspline(xfine[i])    
#        sfint[i] = sfnc(xfine[i])          
        
    plot(xfine,ffine,xfine,fint)    
#    plot(xfine,ffine,xfine,sfint) 
#    print(sfnc)
# test routine        
if __name__ == "__main__":  
    test()   