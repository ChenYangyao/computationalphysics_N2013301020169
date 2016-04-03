# Program: this program solves for the projectile motion
#          with air resistance
# it is according to the problem 2.9 in the textbook
# the form of density change of air is according to 'adiabetic model'
# Author :Chen Yangyao       Last Modify: 20160403

import math     # import packages
import pylab as pl
import numpy as np

# class PROJECTILE solves for projectile motion with air resistance
# where : 
#             x0, y0 : the connon's coordinates
#             v0, theta0 : shell's initial velocity and angle
#             B_m, a, alpha, T0 : coefficient of air resistance
#             g, dt: gravity acceleration and compute step size
class PROJECTILE(object):     
    def __init__(self,_x0,_y0,_v0,_theta0,_B_m,_a,_alpha,_T0,_g,_dt):
        self.x= [float(_x0)]
        self.y= [float(_y0)]
        self.vx= _v0*math.cos(_theta0)
        self.vy= _v0*math.sin(_theta0)
        self.v0, self.theta0= _v0, _theta0
        self.v= _v0
        self.B= float(_B_m)        
        self.a, self.alpha, self.T0= float(_a), float(_alpha), float(_T0)
        self.g=_g
        self.dt= float(_dt)
    def calculate(self):
        while True :
            if self.y[-1]<-1E-4 :
                self.yi= self.y[-1]+self.vy*self.dt
                self.xi= self.x[-1]+self.vx*self.dt
                self.gama= -self.y[-1]/self.yi
                self.x.append((self.x[-1]+self.gama*self.xi)/(self.gama+1.))
                self.y.append(0.)
                break
            self.x.append(self.vx*self.dt+self.x[-1])
            self.y.append(self.vy*self.dt+self.y[-1])
            self.vx= -self.B*((1-self.a*self.y[-2]/self.T0)**self.alpha)\
                           *self.vx*self.v*self.dt+self.vx
            self.vy= -self.B*((1-self.a*self.y[-2]/self.T0)**self.alpha)\
                           *self.vy*self.v*self.dt+self.vy-self.g*self.dt
            self.v= math.sqrt(self.vx**2+self.vy**2)
    def store_max(self,x,v,theta):          # x,v,theta: to store the max range
        x.append(self.x[-1])                #            of the projectile
        v.append(self.v0)
        theta.append(self.theta0*180./np.pi)
    def plot(self):                         # plot the orbit
        pl.plot(self.x,self.y,'--')
        pl.annotate(r'$\theta$'+'= '+str(self.theta0*180./np.pi)+r'$^{o}$',
            xy=(self.x[int(len(self.x)/5)],self.y[int(len(self.x)/5)]), 
            fontsize=12)
class PROJECTILE_DT(PROJECTILE):
    def plot(self):
        pl.plot(self.x,self.y,'--',label='dt= %.2f s'%self.dt)
        
# class PROJECTILE_NONE solves for projectile motion WITHOUT air resistance
# where also: 
#             x0, y0 : the connon's coordinates
#             v0, theta0 : shell's initial velocity and angle
#             g, dt: gravity acceleration and compute step size                    
class PROJECTILE_NONE(PROJECTILE):
    def calculate(self):
        while True :
            if self.y[-1]<-1E-4 :
                self.yi= self.y[-1]+self.vy*self.dt
                self.xi= self.x[-1]+self.vx*self.dt
                self.gama= -self.y[-1]/self.yi
                self.x.append((self.x[-1]+self.gama*self.xi)/(self.gama+1.))
                self.y.append(0.)
                break
            self.x.append(self.vx*self.dt+self.x[-1])
            self.y.append(self.vy*self.dt+self.y[-1])
            self.vx= self.vx
            self.vy= self.vy-self.g*self.dt
    def plot(self):
        pl.plot(self.x,self.y,'--')
        pl.annotate(r'$\theta$'+'= '+str(self.theta0*180./np.pi)+r'$^{o}$',
            xy=(self.x[int(len(self.x)*1/2)],self.y[int(len(self.x)*1/2)]), 
            fontsize=12)
            
# class FIND_MAX give the max range of the projectile 
# where:
#            x: storage of range according to different v0 and theta0
#            v: different v0
#            theta: different theta
class FIND_MAX(object):
    def __init__(self,_x,_v,_theta):
        self.x= _x[:]
        self.v= _v[:]
        self.theta= _theta[:]
        self.xmax= 0
        self.xmax_posi= 0
    def calculate(self):
        for i in range(len(self.x)):
            if self.x[i]>self.xmax :
                self.xmax = self.x[i]
                self.xmax_posi = i

# initialize var. & para.
# initialize var. & para.
# where : 
#             x0, y0 : the connon's coordinates
#             v0, theta0 : shell's initial velocity and angle
#             B_m, a, alpha, T0 : coefficient of air resistance
#             g, dt: gravity acceleration, compute step size
x0, y0= 0., 0.
v0= 700.
B_m, a, alpha, T0, g= 4E-5, 6.5E-3, 2.5, 25.+273.5, 9.8
dt= 0.5
figure= pl.figure(figsize=(18,5))

# plot the motion with air resistance
# where:    theta0= pi/4 rad,   v0= 700m/s,   dt= 0.5s 
pl.subplot(131)
pl.xlim(0,30000)
pl.ylim(0,15000)
pl.xlabel('x(m)',fontsize=15)
pl.ylabel('y(m)',fontsize=15)
comp= PROJECTILE(x0, y0, v0, np.pi/4, B_m, a, alpha, T0, g, 0.5)
comp.calculate()
comp.plot()
pl.annotate(r'$v_{0}$'+'= '+'700.00 m/s',
            xy=(1000,13000), fontsize=12
            )            

# plot the motion with diffrent dt
# where:    theta0= pi/4 rad,   v0= 700m/s,   dt= 20, 10, 5, 0.5, 0.1s
pl.subplot(132)
pl.xlim(0,30000)
pl.ylim(0,15000)
pl.xlabel('x(m)',fontsize=15)
for dt_range in [20.,10.,5.,0.5,0.1] :
    comp= PROJECTILE_DT(x0, y0, v0, np.pi/4, B_m, a, alpha, T0, g, dt_range)
    comp.calculate()
    comp.plot()
pl.annotate(r'$v_{0}$'+'= '+'700.00 m/s',
            xy=(1000,13000), fontsize=12
            )            
pl.annotate(r'$\theta_{0}= 45$'+r'$^{o}$',
            xy=(1000,12000), fontsize=12
            )
pl.annotate('with different "dt"',
            xy=(1000,11000), fontsize=12
            )
pl.legend(loc='best',fontsize=12)

# plot the motion to compare it with case where there isn't air resistance
# where:    theta0= pi/4 rad,   v0= 700m/s,   dt= 0.5s
pl.subplot(133)
x,v,theta= [], [], []
for thetarange in range(30,70,1) :
    thetarange_rad= thetarange*np.pi/180
    comp= PROJECTILE(x0, y0, v0, thetarange_rad, B_m, a, alpha, T0, g, dt)
    comp.calculate()
    comp.store_max(x,v,theta)
    if thetarange in range(35,75,15):
        comp.plot()
findmax= FIND_MAX(x,v,theta)
findmax.calculate()
print '     max range with resistance     '          # give max range with air resistance
print x[findmax.xmax_posi],'m     ',v[findmax.xmax_posi],'m/s    ',theta[findmax.xmax_posi],'degree'

x,v,theta= [], [], []
for thetarange in range(30,70,1) :
    thetarange_rad= thetarange*np.pi/180
    comp= PROJECTILE_NONE(x0, y0, v0, thetarange_rad, B_m, a, alpha, T0, g, dt)
    comp.calculate()
    comp.store_max(x,v,theta)
    if thetarange in range(35,75,15):
        comp.plot()
findmax= FIND_MAX(x,v,theta)
findmax.calculate()
print '     max range without resistance    '    # give max range WITHOUT air resistance
print x[findmax.xmax_posi],'m     ',v[findmax.xmax_posi],'m/s    ',theta[findmax.xmax_posi],'degree'
pl.xlabel('x(m)',fontsize=15)
pl.xlim(0,60000)
pl.ylim(0,26000)
pl.annotate(r'$v_{0}$'+'= '+'700.00 m/s',
            xy=(40000,22000), fontsize=12
            )            
pl.annotate('with & without \n air resistance',
            xy=(40000,18000), fontsize=12
            )
pl.show(figure)
    

            
    
        
    
    
    
    
    
    