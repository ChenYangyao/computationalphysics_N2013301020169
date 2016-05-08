'''
PROGRAM solar gr bending 
this program solves for the bending of light in GR theory
author: ChenYangyao          Last Modify: 20160508
'''

from numpy import *
import matplotlib.pyplot as plt

class BENDING(object):
    ''' 
    class BENDING solves for bending of light
    where:  
        E: Energy  L: angular momentum  M: mass of central star   
        r0, ph0: initial distance and angle between star and source
        dt, time :time step size and total time
    '''
    def __init__(self, _E=0.22, _L=1., _M=1., _r0=2.5,_ph0=0., _dt=0.0001, _time =20.):
        self.E, self.L, self.M=_E, _L, _M
        self.r, self.ph, self.t=[_r0], [_ph0], [0.]
        self.dt= _dt
        self.time=_time
        self.n= int(_time/_dt)
    def cal(self):
        for i in range(self.n):
            self.t.append(self.t[-1]+self.dt)
            self.r.append((sqrt(self.E**2-(1-2*self.M/self.r[-1])*(self.L**2)/self.r[-1]**2)/self.E*(1-self.M*2/self.r[-1]))*self.dt+self.r[-1])
            self.ph.append((self.L/(self.r[-2]**2)/self.E*(1-2*self.M/self.r[-2]))*self.dt+self.ph[-1])
    def plot_trajectory(self,_ax):
        self.x=array(self.r)*cos(array(self.ph))
        self.y=array(self.r)*sin(array(self.ph))
        _ax.plot(self.x,self.y,'-y',label='E=%.2f'%self.E)
        _ax.plot([self.x[0]],[self.y[0]],'ob',markersize=8)
        _ax.plot([self.x[-1]],[self.y[-1]],'or',markersize=8)
    def find_min(self,_search):
        for i in range(len(self.r)-2):
            if self.r[i]>self.r[i+1] and self.r[i+1]<self.r[i+2]:
                return (i+1)
        

fig=plt.figure(figsize=(10,5))
ax1=plt.subplot(111)
cmp=BENDING()
cmp.cal()
cmp.plot_trajectory(ax1)

mfile=open(r'd:\bending.txt','w')
for i in range(len(cmp.r)):
    print >> mfile,cmp.r[i],cmp.ph[i]
mfile.close() 

plt.show()
        
'''        
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
'''



    
            
            
            