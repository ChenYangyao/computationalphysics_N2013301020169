import numpy as np
import math

import mpl_toolkits.mplot3d 
import matplotlib.pyplot as plt

class BASEBALL(object):
    def __init__(self, _vx0, _vy0, _vz0, _dt= 0.1, _omgx=0,_omgy=0,_omgz=0):
        self.vx, self.vy, self.vz= _vx0, _vy0, _vz0 
        self.v = math.sqrt(_vx0**2+ _vy0**2+ _vz0**2)
        self.B2= 0.0039+ 0.0058/(1.+math.exp((self.v-35)/5))
        self.S0= 4.1E-4
        self.g= 9.8
        self.dt= _dt 
        self.x, self.y, self.z= [0], [1.8], [0]
        self.omgx, self.omgy, self.omgz= _omgx, _omgy, _omgz 
    def calculate(self):
        while True: 
            self.x.append(self.vx*self.dt+self.x[-1])
            self.y.append(self.vy*self.dt+self.y[-1])
            self.z.append(self.vz*self.dt+self.z[-1])
            self.vx, self.vy, self.vz = \
                (-self.B2*self.v*self.vx+ self.S0*self.vy*self.omgz)*self.dt+ self.vx, \
                (-self.g- self.B2*self.v*self.vy+ self.S0*self.vz*self.omgx)*self.dt+ self.vy,\
                (self.S0*self.vx*self.omgy)*self.dt+ self.vz
            self.v= math.sqrt(self.vx**2+self.vy**2+self.vz**2)
            self.B2= 0.0039+ 0.0058/(1.+math.exp((self.v-35)/5))
            if self.y[-1]< 0: 
                break
    def graphics(self,_gra, _omgy):
        _gra.plot(self.z, self.x, self.y, label=r'$\omega _y$ = %.2f rad/s'%_omgy)
        _gra.scatter([self.z[0],self.z[-1]],[self.x[0],self.x[-1]],[self.y[0],self.y[-1]],s=30)
        _gra.text(self.z[-1], self.x[-1]-80, self.y[-1], r'$\omega _y$ = %.2f rad/s'%_omgy,fontsize=10)


fig= plt.figure(figsize=(6,6))
ax = plt.subplot(1,1,1,projection='3d')
for omgy in [-400.,200.,0.,200.,400.]:
    comp= BASEBALL(110*0.4470*math.cos(np.pi/4), 110*0.4470*math.sin(np.pi/4), 0., 0.1,0, omgy, 0)
    comp.calculate()
    comp.graphics(ax, omgy)
ax.set_xlabel('z (m)', fontsize=18)
ax.set_ylabel('x (m)', fontsize=18)
ax.set_zlabel('y (m)', fontsize=18)
ax.set_title('ball with horizontal spin', fontsize=18)
ax.set_xlim(-100,100)
ax.set_ylim(0,180)
ax.set_zlim(0,100)
ax.text(0,0,80,'with air drag', fontsize= 18)
plt.show(fig)    