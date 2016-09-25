# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 10:27:55 2016

@author: mol
"""

from scipy import empty, zeros, delete, pi, asarray, array
from scipy.linalg import solve
from pylab import plot
import pdb
import numpy as np
import matplotlib.pyplot as plt
from Finite_Differences import diff1, diff2, diff3
#from openpyxl import Workbook # https://openpyxl.readthedocs.io/en/default/

class FEBeam(object):
    def __init__(self):
        self.__sections = [] # Empty section list

    def addSection(self,section):    
        self.__sections.append(section) # Add beam sections to section list
        
        
    def solve(self):
        self.__makeGlobals() # Gessamtsteifigkeitsmatrik und -Lastvektor des Balken
        self.__makeBCs()
        K     = self.__K     # Steifigkeitsmatrix erzeugen und lösen
        L     = self.__L     # Lastvektor erzeugen
        b     = solve(K,L)
        self.__displacment(b)
        
        
    def __makeGlobals(self): 
        elements = []
        for section in self.__sections:
            elements.extend(section.elements)    
        
        ie    = len(elements)
        K     = zeros([4+(ie-1)*2,4+(ie-1)*2])
        L     = zeros([4+(ie-1)*2])
        x     = zeros([2+(ie-1)])
        i     = 0    
        for element in elements:
            K[i:i+4,i:i+4]     = K[i:i+4,i:i+4] + element.K
            L[i:i+4]           = L[i:i+4] + element.L
            x[i/2]             = element.x1
            i                  = i+2
        x[-1] = elements[-1].x2
        
        self.__K     = K
        self.__Kstar = K
        self.__L     = L
        self.__Lstar = L
        self.__x     = x
        self.__globalsize = 4+(ie-1)*2
        
        
    def __makeBCs(self): # Globale Randbedingungen bauen
# führt zum entfernen von Zeilen und Spalten in K und von Zeilen in L
        K = self.__K
        L = self.__L
        bcindex  = []
        bcvalue  = []
        
        rem = 0
        
        for isec,section in enumerate(self.__sections):
            ie = len(section.elements)
            Ksize = 4+(ie-1)*2
            if isec==0:
                bcindex.append(0) 
                bcvalue.append([section.bc1[0],"F"])
                bcindex.append(1) 
                bcvalue.append([section.bc1[1],"M"])
                bcindex.append(Ksize-2)
                bcvalue.append([section.bc2[0],"F"])  
                bcindex.append(Ksize-1)
                bcvalue.append([section.bc2[1],"M"])
                ind = Ksize-1
            else:
                ind = ind + Ksize-2
                bcindex.append(ind-1)
                bcvalue.append([section.bc2[0],"F"]) 
                bcindex.append(ind)
                bcvalue.append([section.bc2[1],"M"])      
            
        for iind,index in enumerate(bcindex): 
            if bcvalue[iind][0] == 0.:
                ind = index - rem
                K = delete(K,ind,0)
                K = delete(K,ind,1)
                L = delete(L,ind,0)
                rem = rem + 1                
            if isinstance(bcvalue[iind][0],tuple):
                ind = index - rem
                if bcvalue[iind][1] == "F":
                    cwF  = bcvalue[iind][0][0]
                    cfiF = bcvalue[iind][0][1]
                    K[ind,ind] = K[ind,ind] + cwF
                    K[ind+1,ind] = K[ind+1,ind] + cfiF
                if bcvalue[iind][1] == "M":
                    cwM  = bcvalue[iind][0][1]
                    cfiM = bcvalue[iind][0][0]
                    K[ind,ind] = K[ind,ind] + cfiM
                    K[ind-1,ind] = K[ind-1,ind] + cwM
                    
        self.__K = K
        self.__L = L                    
        self.__bcindex = bcindex
        self.__bcvalue = bcvalue     
                    
                    
    def __displacment(self,b):
        sol         = zeros(self.__globalsize)
        j = 0
        for i in range(self.__globalsize):
            if (i in self.__bcindex):
                ind = self.__bcindex.index(i)
                if (self.__bcvalue[ind][0] == 0.):
                #sol[i] = 0
                    pass
                else:
                    sol[i] = b[j]
                    j      = j + 1                
            else:
                sol[i] = b[j]
                j      = j + 1 
        self.__w   = zeros(self.__globalsize/2)
        self.__phi = zeros(self.__globalsize/2)
#        self.__M   = zeros(self.__globalsize/2)

        for i in range(int(self.__globalsize/2)):
            self.__w[i]   = sol[2*i] 
            self.__phi[i] = sol[2*i+1] 
            
        self.__solution = sol

    @property                
    def I(self):
        I = empty(int(self.__globalsize/2))
        for i in range(int(self.__globalsize/2)):
            x = self.__x[i]
            for section in self.__sections:
                if section.x1 <= x and x<= section.x2:
                    I[i] = section.amfunc(x - section.x1)
        return I
        

    @property
    def E(self):
        E = empty(int(self.__globalsize/2))
        for i in range(int(self.__globalsize/2)):
            x = self.__x[i]
            for section in self.__sections:
                if section.x1 <= x and x<= section.x2:
                    E[i] = section.E
        return E
        
        
    @property
    def M(self):
        return diff1(self.__phi,self.__x)*self.E*self.I 
        


    @property
    def Q(self):
        return diff2(self.__phi,self.__x)*self.E*self.I 


    @property
    def q(self):
        return diff3(self.__phi,self.__x)*self.E*self.I 
   
            
    @property
    def w(self):
        return self.__w#*self.E*self.I 
      
        
    @property
    def phi(self):
        return self.__phi#*self.E*self.I 
      
        
    @property
    def x(self):
        return self.__x
        
        
    @property
    def Kstar(self):
        return self.__Kstar 
        
        
    @property
    def reactionforces(self):
        P = self.__Kstar.dot(self.__solution)
        R = []
        for i in range(self.__globalsize):
            if (i in self.__bcindex):
                ind = self.__bcindex.index(i)
                if (self.__bcvalue[ind][0] == 0.):
                    R.append(P[i])

        return R        

        

        
class BeamSection(object):
    def __init__(self,x1,x2,BC1,BC2,E=210.*10**9,ny=0.3):
        self.__x1        = x1
        self.__x2        = x2
        self.__length    = x2 - x1
        self.__noelements()
        self.__BC1       = BC1
        self.__BC2       = BC2
        self.__E         = E
        self.__ny        = ny
        self.__G         = self.__E/2./(1. + self.__ny)
        self.__elements  = []
#        self.__make_elements()
#        self.__am_is_set = False
#        self.__as_is_set = False
    
    def make_elements(self):
        dx = self.__length/self.__NoE
        for i in range (self.__NoE):
            x1      = self.__x1  + i*dx
            x2      = self.__x1  + (i+1)*dx
            I       = self.amfunc((i+(i+1))/2.*dx)
            AS      = self.asfunc((x1 + x2)/2. - self.__x1) #self.__AS1 + (self.__AS2 - self.__AS1)/self.__length*((x2 + x1)/2. - self.__x1)
            element = beamelement(x1,x2,self.__E,I,self.__G,AS)
            self.__elements.append(element)
        

    def addLineLoad(self,func):
        for element in self.__elements:
            x1 = element.x1
            x2 = element.x2
            q1 = func(x1)
            q2 = func(x2)
            element.q1 = q1
            element.q2 = q2
            
            
    def addPointLoad(self,F,x0):
        for element in self.__elements:
            x1 = element.x1
            x2 = element.x2
            if x1<x0 and x0<x2:
                element.pointload = (F,x0)
                exit
            if x1==x0:
                if x1 == self.__x1:
                    element.pointload = (F,x0)
                else:
                    element.pointload = (F/2.,x0)
            if x2==x0:
                if x2 == self.__x2:
                    element.pointload = (F,x0)
                else:
                    element.pointload = (F/2.,x0)    
                    
                    
    def addFreeMoment(self,Q,x0):
        for element in self.__elements:
            x1 = element.x1
            x2 = element.x2
            if x1<x0 and x0<x2:
                element.freemoment = (Q,x0)
                exit
            if x1==x0:
                if x1 == self.__x1:
                    element.freemoment = (Q,x0)
                else:
                    element.freemoment = (Q/2.,x0)
            if x2==x0:
                if x2 == self.__x2:
                    element.freemoment = (Q,x0)
                else:
                    element.freemoment = (Q/2.,x0)   
                    
    def addAreaModulus(self,amfunc):
        self.amfunc = amfunc
#        self.__make_elements()
        

    def addShearArea(self,asfunc):
        self.asfunc = asfunc
#        self.__make_elements()     
        
        
    def __noelements(self):
        self.__NoE = max(int(self.__length*10),20)
        
        
    @property
    def noElements(self):
        return self.__NoE    
        
        
    @property
    def L(self):
        return self.__length

    @property
    def x1(self):
        return self.__x1
        
    @property    
    def x2(self):
        return self.__x2    
        

               
    @property                
    def elements(self):
        return self.__elements
                
    @property
    def bc1(self):
        return self.__BC1
        
    @property
    def bc2(self):
        return self.__BC2          
              
    def set_E_Modulus(self,E):
        self.__E = E
        
    def get_E_Modulus(self):
        return self.__E
        
    E = property(get_E_Modulus,set_E_Modulus)

    def set_ny(self,ny):
        self.__ny = ny
        
    def get_ny(self):
        return self.__ny
        
    ny = property(get_ny,set_ny)        
 

       
class beamelement(object):
    def __init__(self,x1,x2,E,I,G,AS):
        self.__l  = (x2 - x1)
        self.__x1 = x1
        self.__x2 = x2
        self.__K  = empty([4,4])
        self.__L  = empty([4])
        self.__E  = E
        self.__I  = I
        self.__G  = G
        self.__AS = AS
        self.__makeK()
        self.__q1 = 0.0
        self.__q2 = 0.0
        self.__f  = 0.0
        self.__m  = 0.0
        self.__xf = (x2 + x1)/2.
        self.__xm = (x2 + x1)/2.
        self.__makeL()
        
    def __set_q1(self,q1):
        self.__q1 = q1
        self.__makeL()
        
    def __get_q1(self):
        return self.__q1

    q1 = property(__get_q1,__set_q1)        

    def __set_q2(self,q2):
        self.__q2 = q2
        self.__makeL()
        
    def __get_q2(self):
        return self.__q2
        
    q2 = property(__get_q2,__set_q2)
    
    def __set_pointload(self,args):
        self.__f  = args[0]
        self.__xf = args[1]
#        print(self.__x1,self.__x2,args[1])
        if not self.__x1 <= args[1] <= self.__x2:
            raise ValueError
        self.__makeL()
        
    def __get_pointload(self,f,x):
        return self.__f,self.__xf        

    pointload = property(__get_pointload,__set_pointload)
    
    def __set_freemoment(self,args):
        self.__m  = args[0]
        self.__xm = args[1]
        if not self.__x1 <= args[1] <= self.__x2:
            raise ValueError
        self.__makeL()
        
    def __get_freemoment(self,m,x):
        return self.__m,self.__xm        

    freemoment = property(__get_freemoment,__set_freemoment)

    def __set_constant_value(self, value):
        raise SyntaxError 
        
    def __makeK(self):
        E    = self.__E
        I    = self.__I 
        G    = self.__G 
        AS   = self.__AS       
        l    = self.__l
        Phiz = 12.*E*I/G/AS/l/l
        Ky   = E*I/l**3/(1 + Phiz)
        self.__K[0,0] = 12.*Ky 
        self.__K[0,1] = -6.*Ky*l
        self.__K[0,2] = -12.*Ky
        self.__K[0,3] = -6.*Ky*l
        self.__K[1,0] = -6.*Ky*l
        self.__K[1,1] = (4. + Phiz)*Ky*l**2 
        self.__K[1,2] = 6.*Ky*l
        self.__K[1,3] = (2. - Phiz)*Ky*l**2
        self.__K[2,0] = -12.*Ky
        self.__K[2,1] = 6.*Ky*l
        self.__K[2,2] = 12.*Ky
        self.__K[2,3] = 6.*Ky*l
        self.__K[3,0] = -6.*Ky*l
        self.__K[3,1] = (2. - Phiz)*Ky*l**2
        self.__K[3,2] = 6.*Ky*l
        self.__K[3,3] = (4. + Phiz)*Ky*l**2 
        
    def __get_K(self):
        return self.__K
        
    K = property(__get_K,__set_constant_value) 
    

    
    def __makeL(self):
        #init
        self.__L[0] = 0.
        self.__L[1] = 0.
        self.__L[2] = 0.
        self.__L[3] = 0.   
        l  = self.__l
        # line load
        q0 = self.__q1
        dq = (self.__q2 - self.__q1)
        f1 = (q0/2. + 3./20.*dq)*l
        m1 = -(q0/12. + 1./30.*dq)*l*l
        f2 = (q0/2. + 7./20.*dq)*l 
        m2 = -(q0/12. + 1./20.*dq)*l*l
        self.__L[0] = f1
        self.__L[1] = m1
        self.__L[2] = f2
        self.__L[3] = m2
        # pointload
        a = self.__xf - self.__x1
        b = self.__x2 - self.__xf
        f  = self.__f        
        f1 = f*(b/l)**2*(1.+2.*a/l)
        m1 = -f*a*(b/l)**2
        f2 = f*(a/l)**2*(1.+2.*b/l)
        m2 = -f*b*(a/l)**2       
        self.__L[0] = self.__L[0] + f1
        self.__L[1] = self.__L[1] + m1 
        self.__L[2] = self.__L[2] + f2
        self.__L[3] = self.__L[3] + m2
        # free moment
        m  = self.__m
        a  = self.__xm - self.__x1  
        b  = self.__x2 - self.__xm
        f1 = -6.*m*a*b/l**3
        f2 = -f1
        m1 = -m*b*(2.*a - b)/l**2
        m2 = -m*a*(2.*b - a)/l**2
#        m2 = m*(3. - 3.*(1. - a/l)**2 - 4.*a/l)
#        f1 = 2.*(m*a/l + m/l)
#        f2 = -f1
#        m1 = -m - m2 - f2*l
        self.__L[0] = self.__L[0] + f1
        self.__L[1] = self.__L[1] + m1
        self.__L[2] = self.__L[2] + f2
        self.__L[3] = self.__L[3] + m2   
        
        
    def __get_L(self):
        return self.__L
        
    L = property(__get_L,__set_constant_value)    


    def __get_x1(self):
        return self.__x1
        
    x1 = property(__get_x1,__set_constant_value)
        
        
    def __get_x2(self):
        return self.__x2    
        
    x2 = property(__get_x2,__set_constant_value)
    
    
    @property                
    def I(self):
        return self.__I


        
class ld_Load(object):
    def __init__(self):
        self.__q0 = -6000.

    def __call__(self,x):
        return  self.__q0
    
    
class AreaModulus(object):
    def __init__(self,r1,r2):
        self.__I = 0.0000167
            
    def __call__(self,x):
        return self.__I
        
        
class AreaShear(object):
    def __init__(self,r1,r2):
        self.__AS = 5./6.*0.388

    def __call__(self,x):
        return self.__AS    



   
if __name__ == "__main__":  
    LB = 4.9
    AS = 5./6.*0.388
    Iy = 0.0000167
    beam = FEBeam()
    BC1 = [0.,1.]
    BC2 = [0.,1.]
    section = BeamSection(0.,LB,BC1,BC2)
    load = ld_Load()
    AS   = AreaShear()
    AM   = AreaModulus()
    section.addAreaModulus(AM)
    section.addShearArea(AS)
    section.make_elements()
    section.addLineLoad(load)
    beam.addSection(section)
    beam.solve()
    w   = beam.w
    x   = beam.x
    phi = beam.phi

    E  = 210*10**9
    wana = zeros(x.size)
    for i in range(x.size):
        wana[i]=-6000.*LB**4/24./E/Iy*(x[i]/LB-2*(x[i]/LB)**3+(x[i]/LB)**4)
    
    plot(x,w,x,wana) 
