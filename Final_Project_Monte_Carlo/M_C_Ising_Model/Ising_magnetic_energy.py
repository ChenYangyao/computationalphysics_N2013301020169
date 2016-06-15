'''
Program : Ising model magnetic energy
        this program solves the Ising model
        it also gives the plot magnetization and energy vs. time
Auther: Chen Yangyao      Last Modify: 20160613
'''
import matplotlib.pylab as plt
from numpy import *
'''
class : ising
    this class solves the problem of 2-d ising model
where:
    H_magnetic: magnetic field
    Temperature: temperature of heat bath
    length: length of 2-d ising cube
'''
class ISING(object):
    def __init__(self, _H_magnetic=0., _Temperature=0.5, _length=20):
        self.length = int(_length)
        self.H = float(_H_magnetic)
        self.T = float(_Temperature)
        self.system = []
        for i in range(self.length):
            self.system.append([1]*self.length)
    def surrounding(self,_posi):    # give the position of surrounding spin in periodic boundary condition
        if 1 <= _posi[0] <= self.length-2 and 1<=_posi[1]<=self.length-2:
            return [_posi[0]-1,_posi[1]],[_posi[0]+1,_posi[1]],[_posi[0],_posi[1]-1],[_posi[0],_posi[1]+1]
        if _posi[0]==0 and (1<=_posi[1]<=self.length-2):
            return [self.length-1,_posi[1]],[_posi[0]+1,_posi[1]],[_posi[0],_posi[1]-1],[_posi[0],_posi[1]+1]
        if _posi[0]==self.length-1 and  (1<=_posi[1]<=self.length-2):
            return [_posi[0]-1,_posi[1]],[0,_posi[1]],[_posi[0],_posi[1]-1],[_posi[0],_posi[1]+1]
        if 1<=_posi[0]<=self.length-2 and _posi[1]==0:
            return [_posi[0]-1,_posi[1]],[_posi[0]+1,_posi[1]],[_posi[0],self.length-1],[_posi[0],_posi[1]+1]
        if 1<=_posi[0]<=self.length-2 and _posi[1]==self.length-1:
            return [_posi[0]-1,_posi[1]],[_posi[0]+1,_posi[1]],[_posi[0],_posi[1]-1],[_posi[0],0]
        if _posi==[0,0]:
            return [self.length-1,0],[1,0],[0,self.length-1],[0,1]
        if _posi==[self.length-1,0]:
            return [self.length-2,0],[0,0],[self.length-1,self.length-1],[self.length-1,1]
        if _posi==[0,self.length-1]:
            return [self.length-1,self.length-1],[1,self.length-1],[0,self.length-2],[0,0]
        if _posi==[self.length-1,self.length-1]:
            return [self.length-2,self.length-1],[0,self.length-1],[self.length-1,self.length-2],[self.length-1,0]           
    def MCstep(self,_system, _H=0., _T=0.5):  # sweep all spins -- Monte Carlo step
        for i in range(self.length):
            for j in range(self.length):
                self.temp1,self.temp2,self.temp3,self.temp4=self.surrounding([i,j])
                self.temp=array([_system[self.temp1[0]][self.temp1[1]],_system[self.temp2[0]][self.temp2[1]],_system[self.temp3[0]][self.temp3[1]],_system[self.temp4[0]][self.temp4[1]]])
                self.delta_e = sum(self.temp)*_system[i][j]*2.+2.*_system[i][j]*_H
                if self.delta_e <=0:
                    _system[i][j]=-_system[i][j]
                else :
                    self.temp=exp(-self.delta_e/_T)
                    self.temp1=random.rand()
                    if self.temp1<=self.temp:
                        _system[i][j]=-_system[i][j]
    def energy_ave(self,_system,_H=0.):
        self.energy=0.
        for i in range(self.length):
            for j in range(self.length):
                self.temp1,self.temp2,self.temp3,self.temp4=self.surrounding([i,j])
                self.temp=array([_system[self.temp1[0]][self.temp1[1]],_system[self.temp2[0]][self.temp2[1]],_system[self.temp3[0]][self.temp3[1]],_system[self.temp4[0]][self.temp4[1]]])
                self.energy=self.energy+(-sum(self.temp)*_system[i][j])/2.- _system[i][j]*_H
        return self.energy/self.length**2
    def plot_system(self,_sys):
        pass

class ENERGY_MAG(object):            # give the plot of energy and magnetization vs temperature
    def __init__(self):
        self.T=list(arange(0.5,2.1,0.1))+list(arange(2.1,2.5,0.01))+list(arange(2.5,5,0.1))
    def calculate(self):
        self.energy=[]
        self.mag=[]
        self.e_var=[]
        self.mag_var=[]
        print self.T
        for i in self.T:
            self.n=5000
            self.temp_e=0
            self.temp_m=0
            self.temp_es=0
            self.temp_ms=0
            self.cmp=ISING(0., i,10)
            for j in range(self.n):
                self.cmp.MCstep(self.cmp.system,0.,i)
                if j > 20 :
                    self.temp_ee=self.cmp.energy_ave(self.cmp.system,0.)
                    self.temp_mm=mean(self.cmp.system)
                    self.temp_e=self.temp_e+self.temp_ee
                    self.temp_m=self.temp_m+self.temp_mm
                    self.temp_es=self.temp_es+self.temp_ee**2
                    self.temp_ms=self.temp_ms+self.temp_mm**2
            self.energy.append(self.temp_e/(self.n-21.))
            self.mag.append(self.temp_m/(self.n-21.))
            self.e_var.append((self.temp_es/(self.n-21.)-(self.temp_e/(self.n-21.))**2)/i**2)
            self.mag_var.append((self.temp_ms/(self.n-21.)-(self.temp_m/(self.n-21.))**2)/i)
        fig=plt.figure(figsize=(12,6))
        ax1=plt.subplot(221)
        ax2=plt.subplot(222)
        ax3=plt.subplot(223)
        ax4=plt.subplot(224)
        ax1.plot(self.T,self.mag,'ob',markersize=5)
        ax2.plot(self.T,self.energy,'or',markersize=5)
        ax3.plot(self.T,self.mag_var,'og',markersize=5)
        ax4.plot(self.T,self.e_var,'ok',markersize=5)
        plt.show()
cmp=ENERGY_MAG()
cmp.calculate()
mfile=open(r"d:\data.txt",'w')
for i in range(len(cmp.T)):
    print >> mfile, cmp.T[i], cmp.energy[i], cmp.mag[i], cmp.e_var[i], cmp.mag_var[i]
mfile.close
                    
                

        
        
            