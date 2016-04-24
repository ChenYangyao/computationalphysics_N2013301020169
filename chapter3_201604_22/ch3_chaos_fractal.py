'''
PROGRAM chaos_fractal
this program give the Poincare section of the chaotic pendulum
also show the fractal structure of this system in Poincare section
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
    def __init__(self,_Fd=1.2, _theta0=0.2, _omgd=2./3., _size=100., _period=400.):
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
            while self.theta[-1]>np.pi:
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
    def plot_phase(self,_ax,_style):         # the phase-space plot
        _ax.plot(self.theta,self.omg,'o',color=_style,markersize=1,label=r'$F_d = $'+' %.1f'%self.Fd)
    def plot_Poincare(self,_ax):         # the Poincare section plot
        self.t_Poincare, self.omg_Poincare, self.theta_Poincare=[],[],[]
        for i in range(int(np.round(self.period))):
            self.t_Poincare.append(self.t[(i+1)*self.size])
            self.omg_Poincare.append(self.omg[(i+1)*self.size])
            self.theta_Poincare.append(self.theta[(i+1)*self.size])
        _ax.plot(self.theta_Poincare,self.omg_Poincare,'or',markersize=2,label=r'$F_d = $'+' %.1f'%self.Fd)
        
# give Poincare section of chaotic pendulum    
fig=plt.figure(figsize=(12,4))
ax1=plt.subplot(131)
plt.xticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi],[r'$-\pi$',r'$-\pi /2$',r'$0$',r'$\pi/2$',r'$\pi$'])
plt.xlim(-3.5,4.)
plt.ylim(-2.,1.)
ax2=plt.subplot(132)
plt.xlim(1.5,3.4)
plt.ylim(-1.5,-0.4)
plt.title('Poincare Section - Fractal: '+r'$\omega$'+'  versus  '+r'$\theta$',fontsize=18)
ax3=plt.subplot(133)
plt.xlim(3.03,3.06)
plt.ylim(-0.86,-0.81)
# ax1,ax2,ax3: different axes range
cal=CHAOS(1.2)
cal.calculate()
cal.plot_Poincare(ax1)
cal=CHAOS(1.2,0.2,2./3.,100.,3000.)
cal.calculate()
cal.plot_Poincare(ax2)
cal=CHAOS(1.2,0.2,2./3.,100.,20000.)
cal.calculate()
cal.plot_Poincare(ax3)

plt.show()
      
        
        
