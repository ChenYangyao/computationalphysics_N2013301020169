'''
PROGRAM chaos_theta_t
this program solve for the chaotic pendulum 
Author: Chen Yangyao      Last Modify: 20160424
'''
# import packages needed
import matplotlib.pyplot as plt
import numpy as np
# class CHAOS solves for the chaotic pendulum
# the equation considers both damping, driving, and nonlinearity
# where:        Fd, omgd- amplitude and frequency of driving force
#               size,period - number of steps in an period of driving force, total computing number of period
#               theta0 - initial angle position
#               omg0 =0 -initial angular velocity will be zero
class CHAOS(object):
    def __init__(self,_Fd=1.2, _theta0=0.2, _omgd=2./3., _size=100., _period=4.):
        self.theta, self.omg, self.t=[_theta0], [0.0], [0.0]
        self.size, self.period= int(round(_size)), int(round(_period))
        self.dt=(2.*np.pi/_omgd)/self.size
        self.time=(2.*np.pi/_omgd)*self.period
        self.n=int(np.round(self.time/self.dt))
        self.g, self.l, self.q=9.8, 9.8, 1./2.
        self.Fd, self.omgd=_Fd, _omgd 
    def calculate(self):       # use fourth-order Runge-Kutta method to solve the chaotic pendulum
        for i in range(self.n):
            self.t1,self.t2,self.t3,self.t4=self.t[-1],self.t[-1]+self.dt/2.,self.t[-1]+self.dt/2.,self.t[-1]+self.dt
            self.omg1=self.omg[-1]
            self.theta1=self.theta[-1]
            self.omg2=self.omg1+(-self.g/self.l*np.sin(self.theta1)-self.q*self.omg1+ \
                            self.Fd*np.sin(self.omgd*self.t1))*self.dt/2.
            self.theta2=self.theta1+self.omg1*self.dt/2.
            self.omg3=self.omg1+(-self.g/self.l*np.sin(self.theta2)-self.q*self.omg2+ \
                            self.Fd*np.sin(self.omgd*self.t2))*self.dt/2.
            self.theta3=self.theta1+self.omg2*self.dt/2.
            self.omg4=self.omg1+(-self.g/self.l*np.sin(self.theta3)-self.q*self.omg3+ \
                            self.Fd*np.sin(self.omgd*self.t3))*self.dt
            self.theta4=self.theta1+self.omg3*self.dt
            self.t.append(self.t4)
            self.omg.append(self.omg1+ \
                            1./6.*(-self.g/self.l)*(np.sin(self.theta1)+2.*np.sin(self.theta2)+2.*np.sin(self.theta3)+np.sin(self.theta4))*self.dt+ \
                            1./6.*(-self.q)*(self.omg1+2.*self.omg2+2.*self.omg3+self.omg4)*self.dt +\
                            1./6.*self.Fd*(np.sin(self.omgd*self.t1)+2.*np.sin(self.omgd*self.t2)+2.*np.sin(self.omgd*self.t3)+np.sin(self.omgd*self.t4))*self.dt
                            )
            self.theta.append(self.theta1+1./6.*(self.omg1+2.*self.omg2+2.*self.omg3+self.omg4)*self.dt)
            while self.theta[-1]>np.pi:     # reset the angle to keep it in range [-pi,pi]
                self.theta[-1]=self.theta[-1]-2.*np.pi
            while self.theta[-1]<-np.pi:
                self.theta[-1]=self.theta[-1]+2.*np.pi
    def calculate_allangle(self):           # calculate, but don't reset the angle to keep it in range [-pi,pi]
        for i in range(self.n):
            self.t1,self.t2,self.t3,self.t4=self.t[-1],self.t[-1]+self.dt/2.,self.t[-1]+self.dt/2.,self.t[-1]+self.dt
            self.omg1=self.omg[-1]
            self.theta1=self.theta[-1]
            self.omg2=self.omg1+(-self.g/self.l*np.sin(self.theta1)-self.q*self.omg1+ \
                            self.Fd*np.sin(self.omgd*self.t1))*self.dt/2.
            self.theta2=self.theta1+self.omg1*self.dt/2.
            self.omg3=self.omg1+(-self.g/self.l*np.sin(self.theta2)-self.q*self.omg2+ \
                            self.Fd*np.sin(self.omgd*self.t2))*self.dt/2.
            self.theta3=self.theta1+self.omg2*self.dt/2.
            self.omg4=self.omg1+(-self.g/self.l*np.sin(self.theta3)-self.q*self.omg3+ \
                            self.Fd*np.sin(self.omgd*self.t3))*self.dt
            self.theta4=self.theta1+self.omg3*self.dt
            self.t.append(self.t4)
            self.omg.append(self.omg1+ \
                            1./6.*(-self.g/self.l)*(np.sin(self.theta1)+2.*np.sin(self.theta2)+2.*np.sin(self.theta3)+np.sin(self.theta4))*self.dt+ \
                            1./6.*(-self.q)*(self.omg1+2.*self.omg2+2.*self.omg3+self.omg4)*self.dt +\
                            1./6.*self.Fd*(np.sin(self.omgd*self.t1)+2.*np.sin(self.omgd*self.t2)+2.*np.sin(self.omgd*self.t3)+np.sin(self.omgd*self.t4))*self.dt
                            )
            self.theta.append(self.theta1+1./6.*(self.omg1+2.*self.omg2+2.*self.omg3+self.omg4)*self.dt)    
    def plot_theta(self,_ax):         # the theta(angle)-t plot
        _ax.plot(self.t,self.theta,'-',label=r'$F_d = $'+' %.2f'%self.Fd)
    def plot_omg(self,_ax):           # the omega(angular velocity)-t plot
        _ax.plot(self.t,self.omg,'-',label=r'$F_d = $'+' %.2f'%self.Fd)
    def plot_phase(self,_ax):         # the phase-space plot
        _ax.plot(self.theta,self.omg,'-',label=r'$F_d = $'+' %.2f'%self.Fd)
    def plot_Poincare(self,_ax):         # the Poincare section plot
        self.t_Poincare, self.omg_Poincare, self.theta_Poincare=[],[],[]
        for i in range(int(np.round(self.period))):
            self.t_Poincare.append(self.t[(i+1)*self.size])
            self.omg_Poincare.append(self.omg[(i+1)*self.size])
            self.theta_Poincare.append(self.theta[(i+1)*self.size])
        _ax.scatter(self.theta_Poincare[300:], self.omg_Poincare[300:], s=8,label=r'$F_d = $'+' %.2f'%self.Fd)

# class CHAOS_VIA solves two identical pendulum systems
#                    but initial angle has a difference 1E-3
# where:        Fd- amplitude of driving force
#               theta01,theta02 - initial angle of two pendulum 
class CHAOS_VIA(object):
    def __init__(self,_Fd=1.2,_theta01=0.2,_theta02=0.2-1E-3):
        self.Fd=_Fd 
        self.theta0=[_theta01,_theta02]
    def calculate(self):
        self.cal=CHAOS(self.Fd, self.theta0[0],2./3.,100,8)
        self.cal.calculate()
        self.t=self.cal.t
        self.theta1=self.cal.theta
        self.cal=CHAOS(self.Fd, self.theta0[1],2./3.,100,8)
        self.cal.calculate()
        self.theta2=self.cal.theta
        self.theta=np.array(self.theta1)-np.array(self.theta2)
        self.theta=np.abs(self.theta)
    def plot_via(self,_ax):
        _ax.semilogy(self.t, self.theta,'-r',label=r'$F_d = $'+' %.2f'%self.Fd)
        
        
# give figures of chaotic pendulum
fig=plt.figure(figsize=(10,8))
ax1=plt.subplot(321)
plt.title(r'$\theta$'+'  versus  '+r'time',fontsize=18)
ax2=plt.subplot(323)
plt.yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],[r'$-\pi$',r'$-\pi /2$',r'$0$',r'$\pi/2$',r'$\pi$'])
ax2.set_ylim(-3.5,3.5)
ax3=plt.subplot(325)
plt.yticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],[r'$-\pi$',r'$-\pi /2$',r'$0$',r'$\pi/2$',r'$\pi$'])
ax3.set_ylim(-3.5,3.5)

# ax1,ax2,ax3 : theta-t plot
cal=CHAOS(0.,0.2,2./3.,100.,4.)  #Fd=0
cal.calculate()
cal.plot_theta(ax1)
cal=CHAOS(0.5,0.2,2./3.,100.,4.) #Fd=0.5
cal.calculate()
cal.plot_theta(ax2)
cal=CHAOS(1.2,0.2,2./3.,100.,4.) #Fd=1.2
cal.calculate()
cal.plot_theta(ax3)

# ax4,ax5,ax6 : theta-difference of two pendulums,with different initial angle
ax4=plt.subplot(322)
plt.title(r'$\Delta \theta$'+'  versus  '+r'time',fontsize=18)
ax5=plt.subplot(324)
ax6=plt.subplot(326)
cal=CHAOS_VIA(0.0)  #Fd=0
cal.calculate()
cal.plot_via(ax4)
cal=CHAOS_VIA(0.5)  #Fd=0.5
cal.calculate()
cal.plot_via(ax5)
cal=CHAOS_VIA(1.2)  #Fd=1.2
cal.calculate()
cal.plot_via(ax6)

plt.show()

        
        
        
        