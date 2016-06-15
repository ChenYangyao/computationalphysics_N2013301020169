'''
Program : Ising model fluctuation
        this program solves the Ising model
        it also gives the plot magnetization and energy vs. time (given temperature)
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
    def plot_spinvstime(self,_times):         # plot spins energy and magnetization vs. time
        self.time=linspace(0,_times,_times+1)
        self.mag=[mean(self.system)]
        self.ene=[self.energy_ave(self.system,self.H)]
        for i in range(_times):
            self.MCstep(self.system, self.H, self.T)
            self.mag.append(mean(self.system))
            self.ene.append(self.energy_ave(self.system,self.H))
        fig=plt.figure(figsize=(13,12))
        ax1=plt.axes([0.1,0.2,0.35,0.7])
        ax2=plt.axes([0.55,0.2,0.35,0.7])
        ax1.plot(self.time,self.mag,'-b')
        ax2.plot(self.time,self.ene,'-r')
        ax1.set_xlabel("Time / M-C step",fontsize=15)
        ax1.set_ylabel("Magnetization",fontsize=15)
        ax2.set_xlabel("Time / M-C step",fontsize=15)
        ax2.set_ylabel("Energy",fontsize=15)
        ax1.set_title("Fluctuation: Magnetization vs. time",fontsize=18)
        ax2.set_title("Fluctuation: Energy vs. time",fontsize=18)
        ax1.set_ylim(-1.5,1.5)
        ax2.set_ylim(-3,3)
        plt.show()
cmp=ISING(0,0.5)
cmp.plot_spinvstime(200)
    
        
        
            