'''
program anharmonic_a=1_or_3
Author: Chen Yangyao   Last Modify:20160415
'''
# basic packages needed
#     matplotlib - for plot 2D lines
#     numpy - for math 
import matplotlib.pyplot as plt
import numpy as np
# class CROMER
# use RUNGE-KUTTA METHOD to solve the pendulum
# para. & var.  :    
#             theta0, omg0, t0 - initial angel & angular velocity & time
#             l, g, dt, time - length of string, gravity acceleration, time step size & total time 
#             the time period of this pendulum - 1(s)   
class CROMER(object):
    def __init__(self, _theta0=10., _alpha=1.,_omg0=0, _t0=0., _dt=0.001, _time=50.):
        self.theta, self.omg, self.t = [_theta0], [_omg0], [_t0]
        self.alpha=_alpha
        self.dt, self.time, self.n= _dt, _time, int(_time/_dt)
    def calculate(self):
        for i in range(self.n):                  #  RUNGE-KUTTA METHOD
            self.t1,self.t2,self.t3,self.t4=self.t[-1],self.t[-1]+self.dt/2.,self.t[-1]+self.dt/2.,self.t[-1]+self.dt
            self.omg1=self.omg[-1]
            self.theta1=self.theta[-1]
            self.omg2=self.omg1-((self.theta1)**(self.alpha))*self.dt/2.
            self.theta2=self.theta1+self.omg1*self.dt/2.
            self.omg3=self.omg1-(self.theta2**(self.alpha))*self.dt/2.
            self.theta3=self.theta1+self.omg2*self.dt/2.
            self.omg4=self.omg1-(self.theta3**(self.alpha))*self.dt
            self.theta4=self.theta1+self.omg3*self.dt
            self.t.append(self.t[-1]+self.dt)
            self.omg.append(self.omg[-1]-1./6.*((self.theta1)**self.alpha+2.*(self.theta2)**self.alpha  \
                                                    +2.*(self.theta3)**self.alpha+(self.theta4)**self.alpha)*self.dt)
            self.theta.append(self.theta[-1]+1./6.*(self.omg1+2.*self.omg2+2.*self.omg3+self.omg4)*self.dt)
            if (i>500 and np.abs(self.theta[-1]-self.theta[0])<  0.001) or \
                    (i>500 and np.abs(self.theta[-1]-self.theta[0])<np.abs(self.theta[0]-self.theta[1]) ):
                break
        return self.t[-1]-self.t[0]
    def calculate_plot(self,_ax,style):
        for i in range(self.n):
            self.t.append(self.t[-1]+self.dt)
            self.omg.append(self.omg[-1]-((self.theta[-1])**(self.alpha))*self.dt)
            self.theta.append(self.theta[-1]+self.omg[-1]*self.dt)
        _ax.plot(self.t, self.theta, '-'+style,label=r'$\alpha =$'+'%d'%self.alpha)
        
# plot the graphics
n=30
theta=np.linspace(0.2,1.,n)
fig=plt.figure(figsize=(11,6))     # t-theta relation alpha=1,3
ax2=plt.subplot(122) 
ax1=plt.subplot(121)
cal=CROMER(1.,1.)
cal.calculate_plot(ax1,'r')
cal=CROMER(1.,3.)
cal.calculate_plot(ax1,'b')

tau=[]                            # A-theta0 relation alpha=1,3
for i in range(len(theta)):
    cal=CROMER(theta[i],3)
    tau.append(cal.calculate())
    print 'alpha=',3,'theta0=',theta[i],'tau=',tau[-1]
ax2.scatter(theta,tau,s=30,c='blue',label=r'$\alpha=3$'+',Numerical')

tau=[]
for i in range(len(theta)):
    cal=CROMER(theta[i],1.)
    tau.append(cal.calculate())
    print 'alpha=',1,'theta0=',theta[i],'tau=',tau[-1]
ax2.scatter(theta,tau,s=30,c='red',label=r'$\alpha=1$'+',Numerical')

ax2.legend(fontsize=12,loc='best')
ax2.set_xlabel('Amplitude (rad)',fontsize=14)
ax2.set_ylabel('Period (s)',fontsize=14)
ax2.set_title('Amplitude-Period: Change  '+r'$\alpha$',fontsize=14)
ax1.legend(fontsize=12,loc='best')
ax1.set_xlabel('Time (s)',fontsize=14)
ax1.set_ylabel('Angle (rad)',fontsize=14)
ax1.set_title('Oscillator: Change  '+r'$\alpha$',fontsize=14)
ax1.set_xlim(0,25)
ax1.set_ylim(-1.8,1.8)

theta = np.linspace(0.15,1.2,100)      # A-theta0 relation alpha=1,3, Analytical
tau=100*[2*np.pi]
ax2.plot(theta,tau,'--r',label=r'$\alpha=1$'+',Analytical')
tau=2.62206*np.sqrt(8)*theta**(-1)
ax2.plot(theta,tau,'--b',label=r'$\alpha=3$'+',Analytical')
plt.show(fig)