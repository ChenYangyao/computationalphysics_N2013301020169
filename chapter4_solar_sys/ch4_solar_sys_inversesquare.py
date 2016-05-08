'''
PROGRAM solar sys inversesquare
this program solves for the non-inverse-square gravity
author: ChenYangyao          Last Modify: 20160508
'''
from numpy import *
import matplotlib.pyplot as plt

class INVERSE(object):
    ''' 
    class BINARY solves for system that not satisfy the inverse-square law
    where:
        -beta: index of the force
        e: ellipticity
        m: mass of central star
        dt, time : time step size and total time 
    '''
    def __init__(self, _beta=2.05, _e=0., _m=4*(pi**2), _dt=0.001, _time=10):
        self.m=_m
        self.e=_e
        self.x, self.y=[1.],[0.]
        self.vx, self.vy=[0],[sqrt(_m)*sqrt((1.-_e)/(1.+_e))]
        self.beta=_beta
        self.dt=_dt
        self.time= _time
        self.n=int(_time/_dt)
        print self.x[-1],self.y[-1],self.vx[-1],self.vy[-1]
    def cal(self):       # use Euler-Cromer Method to calculate the trajectory of stars
        for i in range(self.n):
            self.r=sqrt(self.x[-1]**2+self.y[-1]**2)
            self.vx.append(self.vx[-1]+self.dt*(-self.m*self.x[-1]/self.r**(self.beta+1.)))
            self.vy.append(self.vy[-1]+self.dt*(-self.m*self.y[-1]/self.r**(self.beta+1.)))
            self.x.append(self.x[-1]+self.vx[-1]*self.dt)
            self.y.append(self.y[-1]+self.vy[-1]*self.dt)
    def plot_trajectory(self,_ax,_style):       # plot the trajectory
        _ax.plot(self.x,self.y,'o'+_style,markersize=0.5,label='e= %.2f'%self.e)
        _ax.plot([self.x[-1]],[self.y[-1]],'o'+_style,markersize=8)
        _ax.plot([0],[0],'or',markersize=20)
    def precession_rate(self):  # calculate the precession rate
        self.x_critical=0
        self.y_critical=0
        self.t_critical=0
        for i in range(len(self.x)-2):
            self.r_i=sqrt(self.x[i]**2+self.y[i]**2)
            self.r_i1=sqrt(self.x[i+1]**2+self.y[i+1]**2)
            self.r_i2=sqrt(self.x[i+2]**2+self.y[i+2]**2)
            if self.r_i<self.r_i1 and self.r_i1>self.r_i2:
                self.x_critical=self.x[i+1]
                self.y_critical=self.y[i+1]
                self.t_critical=self.dt*(i+1)
                break
        self.rate = arctan(self.y_critical/self.x_critical)/self.t_critical
        return self.rate
        
# calculate the trajectory of planet with different ellipticity       
fig=plt.figure(figsize=(10,10)) 
ax1=plt.axes([0.1,0.55,0.35,0.35])
ax2=plt.axes([0.6,0.55,0.35,0.35])
ax3=plt.axes([0.1,0.1,0.35,0.35])
ax4=plt.axes([0.6,0.1,0.35,0.35])

ax1.set_xlim(-1.2,1.2)
ax1.set_ylim(-1.2,1.2)
ax2.set_xlim(-1.2,1.2)
ax2.set_ylim(-1.2,1.2)
ax3.set_xlim(-1.2,1.2)
ax3.set_ylim(-1.2,1.2)
ax4.set_xlim(-0,0.7)


ax1.set_xlabel(r'$x$'+' (AU)',fontsize=18)
ax1.set_ylabel(r'$y$'+' (AU)',fontsize=18)
ax2.set_xlabel(r'$x$'+' (AU)',fontsize=18)
ax2.set_ylabel(r'$y$'+' (AU)',fontsize=18)
ax3.set_xlabel(r'$x$'+' (AU)',fontsize=18)
ax3.set_ylabel(r'$y$'+' (AU)',fontsize=18)
ax4.set_xlabel(r'Ellipticity',fontsize=18)
ax4.set_ylabel(r'Precession rate '+r'$ degree/yr $',fontsize=18)
ax1.set_title(r'$\beta=2.05$'+'  '+r'$e=0$',fontsize=18)
ax2.set_title(r'$\beta=2.05$'+'  '+r'$e=0.15$',fontsize=18)
ax3.set_title(r'$\beta=2.05$'+'  '+r'$e=0.45$',fontsize=18)
ax4.set_title(r'Precession rate',fontsize=18)

cmp=INVERSE(2.05,0)
cmp.cal()
cmp.plot_trajectory(ax1,'b')

mfile=open(r'd:\ee0.txt','w')
for i in range(len(cmp.x)):
    print >> mfile,cmp.x[i],cmp.y[i]
mfile.close()

cmp=INVERSE(2.05,0.15)
cmp.cal()
cmp.plot_trajectory(ax2,'g')

mfile=open(r'd:\ee015.txt','w')
for i in range(len(cmp.x)):
    print >> mfile,cmp.x[i],cmp.y[i]
mfile.close()


cmp=INVERSE(2.05,0.45)
cmp.cal()
cmp.plot_trajectory(ax3,'y')

mfile=open(r'd:\ee045.txt','w')
for i in range(len(cmp.x)):
    print >> mfile,cmp.x[i],cmp.y[i]
mfile.close()

# change the ellipticity and get the corresponding precession rate
e=[]
rate=[]
for i in linspace(0.1,0.6,20):
    cmp= INVERSE(2.05,i)
    cmp.cal()
    e.append(i)
    rate.append(180/pi*cmp.precession_rate())
ax4.plot(e,rate,'oy')
ax4.plot(e,rate,'-r')

plt.show()



    
            
            
            