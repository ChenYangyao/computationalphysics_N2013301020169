from numpy import *
import mpl_toolkits.mplot3d
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import time

class WAVE(object):
    def __init__(self, _n=3000, _steps=1000, _c=300, _x0=0.5):
        self.M=int(_steps)
        self.c=_c
        self.dx=1/float(_steps)
        self.dt=self.dx/float(_c)
        self.r=self.c*self.dt/self.dx
        self.x0=_x0      
        self.n=int(_n)
    def cal(self):
        self.y1=[0]*(self.M+1)
        for i in range(1,self.M,1):
            self.y1[i]=0.03*sin(30.*pi*i*self.dx)
        #for i in range(1,int(0.7*self.M),1):
        #    self.y1[i]=i*self.dx*0.05/0.7
        #for i in range(int(0.7*self.M),self.M):
        #    self.y1[i]=0.05/0.3*(1-i*self.dx)
        self.y2=self.y1[:]
        self.y3=[0.]*(self.M+1)
        for j in range(self.n):
            for i in range(1,self.M,1):
                self.y3[i]=2.*(1-self.r**2)*self.y2[i]-self.y1[i]+self.r**2*(self.y2[i+1]+self.y2[i-1])
            self.y3[0]=0.
            self.y3[self.M]=0.
            self.y1=self.y2[:]
            self.y2=self.y3[:]
            if (mod(j,5)==0):
                self.plot(j*self.dt*1000)
    def plot(self,_t):
        fig=plt.figure(figsize=(10,4))
        self.x=linspace(0.,1.,self.M+1)
        self.ax=plt.axes([0.1,0.2,0.8,0.7])
        self.ax.plot(self.x,self.y2,'-',label='time = %.2f (ms)'%_t)
        self.ax.set_xlim(-0.0,1.0)
        self.ax.set_ylim(-0.06,0.08)
        self.ax.set_xlabel("x (m)",fontsize=18)
        self.ax.set_ylabel("y (m)",fontsize=18)
        self.ax.set_title('Wave on a string: Sine',fontsize=18)
        self.ax.legend(loc='best',fontsize=14)
        plt.savefig(r'd:\save'+'%.2f'%_t+r'.jpg')

cmp=WAVE()
cmp.cal()     
        