
from numpy import *
import matplotlib.pyplot as plt

class BINARY(object):
    ''' 
    class BINARY solves for binary star system
    '''
    def __init__(self, _m1=4*pi, _m2=2*pi, _v10=[0,1./3.*sqrt(6.*pi)], _v20=[0,-2./3.*sqrt(6.*pi)], _dt=0.0001, _time=5):
        self.m1=_m1 
        self.m2=_m2
        self.x1, self.y1= [0.5],[0]
        self.vx1, self.vy1= [_v10[0]],[_v10[1]]
        self.x2, self.y2= [-0.5],[0]
        self.vx2, self.vy2=[_v20[0]], [_v20[1]]
        self.dt=_dt
        self.time= _time
        self.n=int(_time/_dt)
        print self.x1[-1],self.y1[-1],self.vx1[-1],self.vy1[-1]
    def cal(self):
        for i in range(self.n):
            self.r=sqrt((self.x1[-1]-self.x2[-1])**2+(self.y1[-1]-self.y2[-1])**2)
            self.vx1.append(self.vx1[-1]+self.dt*(-self.m2*(self.x1[-1]-self.x2[-1])/self.r**3))
            self.vy1.append(self.vy1[-1]+self.dt*(-self.m2*(self.y1[-1]-self.y2[-1])/self.r**3))
            self.vx2.append(self.vx2[-1]+self.dt*(-self.m1*(self.x2[-1]-self.x1[-1])/self.r**3))
            self.vy2.append(self.vy2[-1]+self.dt*(-self.m1*(self.y2[-1]-self.y1[-1])/self.r**3))
            self.x1.append(self.x1[-1]+self.dt*self.vx1[-1])
            self.y1.append(self.y1[-1]+self.dt*self.vy1[-1])
            self.x2.append(self.x2[-1]+self.dt*self.vx2[-1])
            self.y2.append(self.y2[-1]+self.dt*self.vy2[-1])
    def plot_trajectory(self,_ax):
        _ax.plot(self.x1,self.y1,'-b',label='star '+r'$M_1$')
        _ax.plot(self.x2,self.y2,'-r',label='star '+r'$M_2$')
        _ax.plot([self.x1[-1]],[self.y1[-1]],'ob',markersize=8)
        _ax.plot([self.x2[-1]],[self.y2[-1]],'or',markersize=8)
        
fig=plt.figure(figsize=(10,5)) 
ax1=plt.axes([0.1,0.1,0.35,0.7])
ax2=plt.axes([0.6,0.1,0.35,0.7])
ax1.set_xlim(-0.8,1.4)
ax1.set_ylim(-1.2,1.0)
ax2.set_xlim(-1.,1.2)
ax2.set_ylim(-1.4,0.8)
ax1.set_xlabel(r'$x$'+' (AU)',fontsize=18)
ax1.set_ylabel(r'$y$'+' (AU)',fontsize=18)
ax2.set_xlabel(r'$x$'+' (AU)',fontsize=18)
ax2.set_ylabel(r'$y$'+' (AU)',fontsize=18)
ax1.set_title('Circular orbits',fontsize=18)
ax2.set_title('Elliptical orbits',fontsize=18)
ax1.text(-0.7,0.8,'$M_1/M_2=2$',fontsize=18)
ax2.text(-0.9,0.6,'$M_1/M_2=2$',fontsize=18)

cmp=BINARY()
cmp.cal()
cmp.plot_trajectory(ax1)

cmp=BINARY(4*pi,2*pi,[1./3.*sqrt(6*pi)*sin(pi/6.),1./3.*sqrt(6*pi)*cos(pi/6.)],[-2./3.*sqrt(6*pi)*sin(pi/6.),-2./3.*sqrt(6*pi)*cos(pi/6.)])
cmp.cal()
cmp.plot_trajectory(ax2)

ax1.legend(loc='lower right')
ax2.legend(loc='lower right')

mfile=open(r'd:\binary.txt','r+')
for i in range(len(cmp.x1)):
    print >> mfile,cmp.x1[i],cmp.y1[i],cmp.x2[i],cmp.y2[i]
mfile.close()

plt.show()



    
            
            
            