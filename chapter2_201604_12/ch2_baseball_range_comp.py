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
                self.gama= -self.y[-2]/self.y[-1]
                self.x[-1],self.y[-1]= \
                    (self.x[-2]+self.gama*self.x[-1])/(self.gama+1.),0
                break
        return self.x[-1]
    def graphics(self,_gra, _omgz):
        _gra.plot(self.x, self.y, label=r'$\omega _z$ = %.2f rad/s'%_omgz)
        _gra.scatter([self.x[0],self.x[-1]],[self.y[0],self.y[-1]],s=30)

        
        
class BASEBALL_NONFRIC(BASEBALL):   
    def calculate(self):
        while True: 
            self.x.append(self.vx*self.dt+self.x[-1])
            self.y.append(self.vy*self.dt+self.y[-1])
            self.z.append(self.vz*self.dt+self.z[-1])
            self.vx, self.vy, self.vz = \
                self.vx, \
                -self.g*self.dt+ self.vy, \
                self.vz
            if self.y[-1]< 0:
                self.gama= -self.y[-2]/self.y[-1]
                self.x[-1],self.y[-1]= \
                    (self.x[-2]+self.gama*self.x[-1])/(self.gama+1.),0
                break
        return self.x[-1]
    def graphics(self,_gra,_dt):
        _gra.plot(self.x, self.y, '--',label= 'dt = %.2f s'%_dt)
        _gra.scatter([self.x[0],self.x[-1]],[self.y[0],self.y[-1]],s=30)
        
def findmax(_t):
    maxval=0 
    maxposi=0 
    for i in range(len(_t)):
        if _t[i] > maxval: 
            maxval = _t[i]
            maxposi = i 
    return [maxposi,maxval]

theta = np.linspace(np.pi/180, np.pi*89/180, 89)
r_drag = []
r_nondrag = []
for j in theta:
    print 'theta= ', j*180/np.pi
    comp= BASEBALL_NONFRIC(110*0.4470*math.cos(j), 110*0.4470*math.sin(j), 0.)
    r_nondrag.append(comp.calculate())
    comp= BASEBALL(110*0.4470*math.cos(j), 110*0.4470*math.sin(j), 0.)
    r_drag.append(comp.calculate())
    
posi_drag = findmax(r_drag)[0]
posi_nondrag = findmax(r_nondrag)[0]
print 'with drag the max range:'
print theta[posi_drag]*180/np.pi,'deg    ', r_drag[posi_drag],'m   ' 
print 'without drag the max range:'
print theta[posi_nondrag]*180/np.pi,'deg    ', r_nondrag[posi_nondrag],'m   ' 



   

