# -*- coding: utf-8 -*-
"""
Created on: Wed Apr  6 14:58:30 2016
@author: cmsteinbach
mail to: steinbach@tuhh.de
Abstract: Solving the 1D Burgers equation with the semi-descrete approach of KT
in space and Heun RK in time
Source: New High-Resolution Central Schemes for Nonlinear Conservation Laws and
Convection–Diffusion Equations (Kurganov & Tadmor 2000)
Version: 0.1
Changed on: Wed Apr  6 14:58:30 2016
TO DO´s: -
Parameter List Description: -
Public Methods Description: -
"""
from pylab import linspace, zeros, empty, plot
from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
#from mpl_toolkits.mplot3d.axes3d import get_test_data
#import inspect

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
        print ('FV_FD_1D KT:')
        print ('  Python version')
        print ('')
        print ('  Solve the nonlinear advection equation with source term in 1D,')
        print ('  du/dt = -dF(u)/dx + S(u)')
        print ('  over the interval:')
        print ('    0.0 <= x <= 1.0')
        print ('  with fixed boundary conditions, and')
        print ('  with a given initial condition')
        print ('')
        print ('  We use the Kurganov & Tadmor mixed FV/FD method with Huen RK.')

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
        return u*u/2.
        
    def source(self,u):
        my = 100.
        Q  = -my*u*(u-1.)*(u-0.5)
#        print(Q)
        return Q
        
    def timestep(self,u,jt,dt,dx):
        u[0,jt]  = u[0,jt-1]
        u[1,jt]  = u[1,jt-1]
        dx       = self.__dx
        dt       = self.__dt
        ustar    = zeros(self.__nx)
        ustar[0] = u[0,jt-1]
        ustar[1] = u[1,jt-1]
        dudt    = zeros([self.__nx,2])
        for ix in range(2,self.__nx-2):
            src        = self.source(u[ix  ,jt-1])
            rH,lH      = self.fluxcalc(u[ix-2,jt-1],u[ix-1,jt-1],u[ix  ,jt-1],u[ix+1,jt-1],u[ix+2,jt-1])
            dudt[ix,0] = (rH - lH)/dx + src
            ustar[ix]  = u[ix,jt-1] + dudt[ix,0]*dt
        ustar[self.__nx-1] = u[self.__nx-1,jt-1]
        ustar[self.__nx-2] = u[self.__nx-2,jt-1]

        for ix in range(2,self.__nx-2):
            src        = self.source(ustar[ix])
            rH,lH      = self.fluxcalc(ustar[ix-2],ustar[ix-1],ustar[ix],ustar[ix+1],ustar[ix+2])
            dudt[ix,1] = (rH - lH)/dx + src
            u[ix,jt]   = u[ix,jt-1] + 0.5*(dudt[ix,0] + dudt[ix,1])*self.__dt            
            
        u[self.__nx-2,jt] = u[self.__nx-3,jt]
        u[self.__nx-1,jt] = u[self.__nx-2,jt]
        return u    

    def fluxcalc(self,um2,um1,u,up1,up2):       
        rim1 = self.superbee(um2,um1,u)
        ri   = self.superbee(um1,u,up1)
        rip1 = self.superbee(u,up1,up2) 
        ulm  = um1 + rim1/2.*(u   - um1)
        ulp  = u   - ri/2.  *(up1 - u  )
        urm  = u   + ri/2.  *(u   - um1)
        urp  = up1 - rip1/2.*(up1 - u  )
#        print(ulm,ulp,urm,urp)
        lHm  = self.flux(ulm)
        lHp  = self.flux(ulp)            
        rHm  = self.flux(urm)            
        rHp  = self.flux(urp)
        al   = max(abs(self.Jacobian(ulm)),abs(self.Jacobian(ulp)))
        ar   = max(abs(self.Jacobian(urm)),abs(self.Jacobian(urp)))
#        print(al,ar)
        
        lH   = 0.5*(lHm + lHp - al*(ulp - ulm))            
        rH   = 0.5*(rHm + rHp - ar*(urp - urm))
        return lH, rH
        
    def Jacobian(self,u):
        eps = 0.0001
        dF  = self.flux(u + eps)
        dF  = dF - self.flux(u - eps)
        return dF/2./eps
        

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
        return max(0.,min(theta*r,(1.+r)/2.,theta)) 
        
        
    def superbee(self,um1,u,up1):
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
        return max(0.,min(2.*r,1.),min(r,2.))         

            
    def solve(self):
        print ('')
        print ('  Number of nodes NX = %d' % ( self.__nx ))
        print ('  Number of time steps NT = %d' % ( self.__nt ))
        print ('  left BC, right BC = ',self.__uleft, self.__uright)
        for j in range(1,self.__nt):
            self.__u = self.timestep(self.__u,j,self.__dt,self.__dx)
            
        print ('')
        print ('FV_FD_1D_BURGERS_KT_TEST')
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
    
    