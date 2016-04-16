'''
program SIMPLE_pendulum 
Author: Chen Yangyao   Last Modify:20160415
'''
# basic packages needed
#     matplotlib - for plot 2D lines
#     numpy - for math 
import matplotlib.pyplot as plt
import numpy as np

# class EULER
# use EULER METHOD to solve the pendulum
# para. & var.  :    
#             theta0, omg0, t0 - initial angel & angular velocity & time
#             l, g, dt, time - length of string, gravity acceleration, time step size & total time 
#             the time period of this pendulum - 1(s)
class EULER(object):
    def __init__(self, _theta0=10., _omg0=0, _t0=0., _l=9.8/(4*(np.pi)**2),_g=9.8, _dt=0.01, _time=4.):
        self.theta, self.omg, self.t = [_theta0], [_omg0], [_t0]
        self.l, self.g, self.dt, self.time, self.n= _l, _g, _dt, _time, int(_time/_dt)
        self.E = [0.5*((self.l*_omg0)**2+self.g*self.l*(_theta0)**2)]
    def calculate(self):
        for i in range(self.n):
            self.t.append(self.t[-1]+self.dt)
            self.omg.append(self.omg[-1]-self.g/self.l*self.theta[-1]*self.dt)
            self.theta.append(self.theta[-1]+self.omg[-2]*self.dt)
            self.E.append(0.5*((self.l*self.omg[-1])**2+self.g*self.l*(self.theta[-1])**2))
    def plot_theta(self,_ax):
        _ax.plot(self.t, self.theta, '--',label='Euler Method')
    def plot_E(self,_ax):
        _ax.plot(self.t,self.E,'--',label='Euler Method')

# class CROMER
# use EULER-CROMER METHOD to solve the pendulum
# para. & var.  :    
#             theta0, omg0, t0 - initial angel & angular velocity & time
#             l, g, dt, time - length of string, gravity acceleration, time step size & total time 
#             the time period of this pendulum - 1(s)        
class CROMER(EULER):
    def calculate(self):
        for i in range(self.n):
            self.t.append(self.t[-1]+self.dt)
            self.omg.append(self.omg[-1]-self.g/self.l*self.theta[-1]*self.dt)
            self.theta.append(self.theta[-1]+self.omg[-1]*self.dt)
            self.E.append(0.5*((self.l*self.omg[-1])**2+self.g*self.l*(self.theta[-1])**2))
    def plot_theta(self,_ax):
        _ax.plot(self.t, self.theta, '--',label='Euler-Cromer Method')
    def plot_E(self,_ax):
        _ax.plot(self.t,self.E,'--',label='Euler-Cromer Method')
        
# plot :
#        ax1 - time dependence of angel
#        ax2 - time dependence of energy
# both EULER & EULER CROMER METHOD are used
fig= plt.figure(figsize=(10,5))
ax1 = plt.subplot(121)
ax2 = plt.subplot(122)

comp= EULER()
comp.calculate()
comp.plot_theta(ax1)
comp.plot_E(ax2)
comp= CROMER()
comp.calculate()
comp.plot_theta(ax1)
comp.plot_E(ax2)

t=np.linspace(0,4,100)
theta=10*np.cos(2.*np.pi*t)
E= 0.5*9.8*9.8/(4*(np.pi)**2)*(10**2)
ax1.plot(t,theta,'--',label='Analysis')
ax2.plot(t,len(t)*[E],label='Analysis')

ax2.text(0.5,350,'Period:  '+r'$ \tau $'+' = 1s',fontsize=14)

ax1.set_title('Simple Pendulum - Angle',fontsize=14)
ax2.set_title('Simple Pendulum - Energy',fontsize=14)
ax1.set_xlabel('time'+r'$ /\tau $',fontsize=14)
ax1.set_ylabel('Angel (rad)',fontsize=14)
ax2.set_xlabel('time'+r'$ /\tau $',fontsize=14)
ax2.set_ylabel('Energy (J)',fontsize=14)
ax1.legend(fontsize=12,loc='best')
ax2.legend(fontsize=12,loc='best')
plt.show(fig)
        