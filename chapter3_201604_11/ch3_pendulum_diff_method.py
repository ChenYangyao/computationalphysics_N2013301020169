'''
program SIMPLE_pendulum_diff_method
Author: Chen Yangyao   Last Modify:20160415
'''
# basic packages needed
#     matplotlib - for plot 2D lines
#     numpy - for math 
#     time - for counting time
import time
import matplotlib.pyplot as plt
import numpy as np

# class EULER
# use EULER METHOD to solve the pendulum
# para. & var.  :    
#             theta0, omg0, t0 - initial angel & angular velocity & time
#             l, g, dt, time - length of string, gravity acceleration, time step size & total time 
#             the time period of this pendulum - 1(s)
class EULER(object):
    def __init__(self, _theta0=10., _omg0=0., _t0=0., _l=9.8/(4.*(np.pi)**2),_g=9.8, _dt=0.01, _time=20000.):
        self.theta, self.omg, self.t = [_theta0], [_omg0], [_t0]
        self.l, self.g, self.dt, self.time, self.n= _l, _g, _dt, _time, int(_time/_dt)
        self.E = [0.5*((self.l*_omg0)**2+self.g*self.l*(_theta0)**2)]
    def calculate(self):
        for i in range(self.n):
            self.t.append(self.t[-1]+self.dt)
            self.omg.append(self.omg[-1]-self.g/self.l*self.theta[-1]*self.dt)
            self.theta.append(self.theta[-1]+self.omg[-2]*self.dt)
            self.E.append(0.5*((self.l*self.omg[-1])**2+self.g*self.l*(self.theta[-1])**2))
    def plot_E(self,_ax):
        _ax.plot(self.t,self.E,'--',label='Euler')
    def var(self):
        self.diff=np.array(self.E)-self.E[0]
        self.diff=np.sqrt(np.mean(self.diff*self.diff))
        return self.diff
        
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
    def plot_E(self,_ax):
        _ax.plot(self.t,self.E,'--y',label='Euler-Cromer')
        
# class RUNGE_RUNGE_22
# use 2nd Runge-Kutta Method to solve the pendulum
# para. & var.  :    
#             theta0, omg0, t0 - initial angel & angular velocity & time
#             l, g, dt, time - length of string, gravity acceleration, time step size & total time 
#             the time period of this pendulum - 1(s)  
        
class RUNGE_RUNGE_22(EULER):
    def calculate(self):
        for i in range(self.n):
            self.t2=self.t[-1]+self.dt/2.
            self.omg2=self.omg[-1]-self.g/self.l*self.theta[-1]*self.dt/2.
            self.theta2=self.theta[-1]+self.omg[-1]*self.dt/2.
            self.t.append(self.t[-1]+self.dt)
            self.omg.append(self.omg[-1]-self.g/self.l*self.theta2*self.dt)
            self.theta.append(self.theta[-1]+self.omg2*self.dt)
            self.E.append(0.5*((self.l*self.omg[-1])**2+self.g*self.l*(self.theta[-1])**2))
    def plot_E(self,_ax):
        _ax.plot(self.t,self.E,'--',label='2nd-order Runge-Kutta')

# class RUNGE_RUNGE_44
# use 4th Runge-Kutta Method to solve the pendulum
# para. & var.  :    
#             theta0, omg0, t0 - initial angel & angular velocity & time
#             l, g, dt, time - length of string, gravity acceleration, time step size & total time 
#             the time period of this pendulum - 1(s)        
class RUNGE_RUNGE_44(EULER):
    def calculate(self):
        for i in range(self.n):
            self.t1,self.t2,self.t3,self.t4=self.t[-1],self.t[-1]+self.dt/2.,self.t[-1]+self.dt/2.,self.t[-1]+self.dt
            self.omg1=self.omg[-1]
            self.theta1=self.theta[-1]
            self.omg2=self.omg1-self.g/self.l*self.theta1*self.dt/2.
            self.theta2=self.theta1+self.omg1*self.dt/2.
            self.omg3=self.omg1-self.g/self.l*self.theta2*self.dt/2.
            self.theta3=self.theta1+self.omg2*self.dt/2.
            self.omg4=self.omg1-self.g/self.l*self.theta3*self.dt
            self.theta4=self.theta1+self.omg3*self.dt
            self.t.append(self.t[-1]+self.dt)
            self.omg.append(self.omg[-1]+1./6.*(-self.g/self.l)*(self.theta1+2.*self.theta2+2.*self.theta3+self.theta4)*self.dt)
            self.theta.append(self.theta[-1]+1./6.*(self.omg1+2.*self.omg2+2.*self.omg3+self.omg4)*self.dt)
            self.E.append(0.5*((self.l*self.omg[-1])**2+self.g*self.l*(self.theta[-1])**2))
    def plot_E(self,_ax):
        _ax.plot(self.t,self.E,'--',label='4th-order Runge-Kutta')
        
        
# class VERLET
# use Verlet Method to solve the pendulum
# para. & var.  :    
#             theta0, omg0, t0 - initial angel & angular velocity & time
#             l, g, dt, time - length of string, gravity acceleration, time step size & total time 
#             the time period of this pendulum - 1(s)
class VERLET(EULER):
    def calculate(self):
        self.t1,self.t2,self.t3,self.t4=self.t[-1],self.t[-1]+self.dt/2.,self.t[-1]+self.dt/2.,self.t[-1]+self.dt
        self.omg1=self.omg[-1]
        self.theta1=self.theta[-1]
        self.omg2=self.omg1-self.g/self.l*self.theta1*self.dt/2.
        self.theta2=self.theta1+self.omg1*self.dt/2.
        self.omg3=self.omg1-self.g/self.l*self.theta2*self.dt/2.
        self.theta3=self.theta1+self.omg2*self.dt/2.
        self.omg4=self.omg1-self.g/self.l*self.theta3*self.dt
        self.theta4=self.theta1+self.omg3*self.dt
        self.t.append(self.t[-1]+self.dt)
        self.theta.append(self.theta[-1]+1./6.*(self.omg1+2.*self.omg2+2.*self.omg3+self.omg4)*self.dt)
        for i in range(self.n-1):
            self.t.append(self.t[-1]+self.dt)
            self.theta.append(2.*self.theta[-1]-self.theta[-2]-self.g/self.l*self.theta[-1]*(self.dt)**2)
            self.omg.append((self.theta[-1]-self.theta[-3])/(2.*self.dt))
            self.E.append(0.5*((self.l*self.omg[-1])**2+self.g*self.l*(self.theta[-2])**2))
        self.omg.append(self.omg[-1])
        self.E.append(self.E[-1])
    def plot_E(self,_ax):
        _ax.plot(self.t,self.E,'--',label='Verlet')
# function main: use a method, plot the energy line, and give the accuracy and time cost         
def main(_method,_ax):
    start=time.clock()
    cal=_method()
    cal.calculate()
    end=time.clock()
    cal.plot_E(_ax)
    return [1./cal.var(), np.abs(end-start)]

# plot the graphics
#      period number : 20000
n=[1.,2.,3.]
e_diff = []
cp_time = []
fig= plt.figure(figsize=(12,5))
ax1=plt.subplot(121)
ax2=plt.subplot(122)
ax3=ax2.twinx()

temp = main(CROMER,ax1)   # CROMER METHOD
e_diff.append(temp[0])
cp_time.append(temp[1])

temp = main(VERLET,ax1)    #VERLET METHOD
e_diff.append(temp[0])
cp_time.append(temp[1])

temp = main(RUNGE_RUNGE_44,ax1)   # 4-th RUNGE-KUTTA METHOD
e_diff.append(temp[0])
cp_time.append(temp[1])



ax2.scatter(n, e_diff,c='blue',s=80,label=' Accuracy')
ax3.scatter(n,cp_time,c='red',s=80,label='Time Cost')

ax1.set_xlabel('time'+r'$/\tau$',fontsize=14)
ax1.set_ylabel('Energy (J)',fontsize=14)
ax1.set_title('Energy: Using Different Method',fontsize=14)
ax2.set_xlabel('Method No.',fontsize=14)
ax2.set_ylabel('Accuracy = 1/Standard Error (1/J)',fontsize=14)
ax2.set_title('Comparison of Different Method',fontsize=14)
ax3.set_ylabel('Time Cost (s)',fontsize=14)

ax1.legend(loc='best',fontsize=12)
ax2.legend(loc='lower right')
ax3.legend(loc='upper left')
print 1./np.array(e_diff)
print cp_time
plt.show()
        
        
        
        