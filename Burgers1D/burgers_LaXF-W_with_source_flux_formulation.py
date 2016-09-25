# -*- coding: utf-8 -*-
"""
Created on: Wed Apr  6 14:58:30 2016
@author: cst
mail to: your email adress here
Abstract: Short Description here
Description: Detailed Description here
Source: Sources used for this 
Version: 0.1
Changed on: Wed Apr  6 14:58:30 2016
TO DO´s: -
Parameter List Description: -
Public Methods Description: -
"""
from pylab import empty, arange, linspace, tanh, zeros, empty,sin,pi,plot,sign
from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d.axes3d import get_test_data
import inspect

class Burgers(object):
    def __init__(self,dx,dt,nx,nt):
        self.__dx = dx
        self.__dt = dt        
        self.__nx = nx
        self.__nt = nt
        self.__u  = empty([nx,nt])
        self.__fi = zeros([nx,nt])
        self.__cfl = 0.
        print ('')
        print ('FD_1D_LAX_WENDROFF:')
        print ('Python version')
        print ('')
        print ('Solve the Burrgers equation with source term in 1D,')
        print ('du/dt = -dF(u)/dx + S(u)')
        print ('over the interval:')
        print ('0.0 <= x <= 1.0')
        print ('with fixed boundary conditions, and')
        print ('with a given initial condition')
        print ('')
        print ('Numerics: Lax-Wendroff method with flux blended with Lax-Friedrichs')

    def leftbc(self,uleft):
        self.__uleft = uleft

    def rightbc(self,uright):
        self.__uright = uright
        
    def xjump(self,xj):
        self.__xj = xj
        
    def initial(self):
        for i in range(self.__nx):
            if i*self.__dx < self.__xj:
                self.__u[i,0] = self.__uleft
            else:
                self.__u[i,0] = self.__uright
                
    def flux(self,u):
        return u#*u/2.
        
    def source(self,u):
        my = 1500.
        Q  = -my*u*(u-1.)*(u-0.5)
#        print(Q)
        return Q
        
    def timestep(self,u,jt,dt,dx):
        eps     = 0.0000001
        u[0,jt] = u[0,jt-1]
        dx      = self.__dx
        dt      = self.__dt
        for ix in range(1,self.__nx-2):
            lflx     = self.flux(u[ix-1,jt-1])
            cflx     = self.flux(u[ix,jt-1])
            rflx     = self.flux(u[ix+1,jt-1])
            lu       = u[ix-1,jt-1]
            cu       = u[ix  ,jt-1]
            ru       = u[ix+1,jt-1]
            lJF      = (cflx - lflx)/(u[ix  ,jt-1] - u[ix-1,jt-1] + eps)  
            rJF      = (rflx - cflx)/(u[ix+1,jt-1] - u[ix  ,jt-1] + eps)
            lsrc     = self.source(u[ix-1,jt-1])
            csrc     = self.source(u[ix  ,jt-1])
            rsrc     = self.source(u[ix+1,jt-1])            
            JS       = (rsrc - lsrc)/(u[ix+1,jt-1] - u[ix-1,jt-1] + eps)
            src1     = -csrc*(rJF - lJF)/dx - (lJF + rJF)*(rsrc - lsrc)/dx/4.
            src2     = csrc*JS - JS*(rflx - lflx)/dx/2. 
            cfi      = self.minmod(u[ix-1,jt-1],u[ix  ,jt-1],u[ix+1,jt-1])
    
            lH       = lflx - (1.-cfi)*dx/dt*(cu - lu) - (cfi)*lJF*(cflx - lflx)/dx*dt 
            rH       = rflx - (1.-cfi)*dx/dt*(ru - cu) - (cfi)*rJF*(rflx - cflx)/dx*dt
            
            u[ix,jt] = cu - (rH - lH)/2./dx*dt
            u[ix,jt] = u[ix,jt] + csrc*dt + (cfi)*(src1 + src2)*dt**2/2.
            
        u[self.__nx-2,jt] = u[self.__nx-3,jt]
        u[self.__nx-1,jt] = u[self.__nx-2,jt]
        return u               

    def minmod(self,um1,u,up1):
        # calculate r_{i} = \frac{u_i-u_{i-1}}{u_{i+1}-u_i}
        nom   = (u   - um1)
        denom = (up1 - u  )
          # make sure division by 0 does not happen
        if(abs(nom) < 1e-14): # nom = 0
            nom   = 0.
            denom = 1.
        elif (nom > 1e-14 and abs(denom) < 1e-14): # nom > 0 => r = \inf
            nom   = 1e14
            denom = 1.
        elif (nom < -1e-14 and abs(denom) < 1e-14): # nom < 0 => r = 0
            nom   = -1e14
            denom = 1.
        r = nom/denom
        theta = 2.
        return max(0.8,min(theta*r,1.))

    def minmod2(self,um1,u,up1):
        # calculate r_{i} = \frac{u_i-u_{i-1}}{u_{i+1}-u_i}
        nom   = (u   - um1)
        denom = (up1 - u  )
        if(abs(nom) < 1e-14): # nom = 0
            nom   = 0.
            denom = 1.
        elif (nom > 1e-14 and abs(denom) < 1e-14): # nom > 0 => r = \inf
            nom   = 1e14
            denom = 1.
        elif (nom < -1e-14 and abs(denom) < 1e-14): # nom < 0 => r = 0
            nom   = -1e14
            denom = 1.
        r = nom/denom
        theta = 2.
        beta  =  .0
        return max(beta,min(theta*r,1.))

        
    def mc(self,um1,u,up1):
        # calculate r_{i} = \frac{u_i-u_{i-1}}{u_{i+1}-u_i}
        nom   = (u   - um1)
        denom = (up1 - u  )
          # make sure division by 0 does not happen
        if(abs(nom) < 1e-14): # nom = 0
            nom   = 0.
            denom = 1.
        elif (nom > 1e-14 and abs(denom) < 1e-14): # nom > 0 => r = \inf
            nom   = 1e14
            denom = 1.
        elif (nom < -1e-14 and abs(denom) < 1e-14): # nom < 0 => r = 0
            nom   = -1e14
            denom = 1.
        r = nom/denom
        theta = 2.
        return max(0.,min(theta*r,(1.+r)/2.,1.))      
            
    def solve(self):
        print ('')
        print ('  Number of nodes NX = %d' % ( self.__nx ))
        print ('  Number of time steps NT = %d' % ( self.__nt ))
        print ('  left BC, right BC = ',self.__uleft, self.__uright)
        for j in range(1,self.__nt):
            self.__u = self.timestep(self.__u,j,self.__dt,self.__dx)
            
        print ('')
        print ('FD_1D_ADVECTION_LAX_WENDROFF_TEST')
        print ('  Normal end of execution.')
        print ('')


    def plot3D(self):
        X = zeros ( [ self.__nx, self.__nt ] )
        Y = zeros ( [ self.__nx, self.__nt ] )
        Z = zeros ( [ self.__nx, self.__nt ] )
        x = linspace(0.,self.__nx*self.__dx,self.__nx)
        for i in range(self.__nt):
            X[:,i] = x
            Y[:,i] = i*self.__dt
            
        Z = self.__u 
        fig = plt.figure ( )
        ax = fig.gca ( projection = '3d' )
        surf = ax.plot_surface ( X, Y, Z, cmap = cm.coolwarm, \
        linewidth = 0, antialiased = False )
        ax.set_xlabel ( '<--X-->' )
        ax.set_ylabel ( '<--T-->' )
        ax.set_zlabel ( '<--U(X,T)-->' )
        fig.colorbar ( surf, shrink = 0.5, aspect = 10 )
        plt.show ( )          

    def plot2D(self,jt):
        x = linspace(0.,self.__nx*self.__dx,self.__nx)
        plot(x,self.__u[:,jt])
        
    def cfl(self):
        print ('max CFL condition: dt = (%g) <= (%g) = dx / J / dt', self.__cfl)
    

if __name__ == "__main__":
    burgers = Burgers(0.005,0.001,201,500)
    burgers.leftbc(1.)
    burgers.rightbc(0.0)
    burgers.xjump(.3)
    burgers.initial()
    burgers.solve()
    burgers.plot2D(499)
    burgers.plot3D()    
    
    