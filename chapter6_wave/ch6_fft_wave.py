from numpy import *
import mpl_toolkits.mplot3d
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import time

class WAVE(object):
    def __init__(self,_x0=0.5, _n=5000, _steps=600, _c=300):
        self.M=int(_steps)
        self.c=_c
        self.dx=1/float(_steps)
        self.dt=self.dx/float(_c)
        self.r=self.c*self.dt/self.dx
        self.x0=_x0      
        self.n=int(_n)
    def cal_spec(self,_spec=0.05):
        self.spec_x=int(_spec*self.M)
        self.spec_t=[]
        self.spec_y=[]
        self.y1=[0]*(self.M+1)
        for i in range(1,self.M,1):
            self.y1[i]=0.05*exp(-(i*self.dx-0.4)**2/0.15**2)*sin(20.*pi*i*self.dx)
        self.y2=self.y1[:]
        self.y3=[0.]*(self.M+1)
        for j in range(self.n):
            for i in range(1,self.M,1):
                self.y3[i]=2.*(1-self.r**2)*self.y2[i]-self.y1[i]+self.r**2*(self.y2[i+1]+self.y2[i-1])
            self.y3[0]=0.
            self.y3[self.M]=0.
            self.y1=self.y2[:]
            self.y2=self.y3[:]
            self.spec_t.append(j*self.dt)
            self.spec_y.append(self.y2[self.spec_x])
        self.sp=fft.fft(self.spec_y,norm='ortho')
        self.freq=fft.fftfreq(self.n)/self.dt    
    def plot_spec(self):
        fig=plt.figure(figsize=(11,6))
        ax1=plt.axes([0.1,0.1,0.35,0.8])
        ax2=plt.axes([0.55,0.1,0.35,0.8])
        ax1.plot(self.spec_t,self.spec_y,'-',lw=1,label="observer at %.2f"%(self.spec_x/float(self.M)))
        ax2.plot(self.freq, abs(self.sp)**2, '-r', label='power spectrum')
        ax1.legend(loc='best')
        ax2.set_xlim(0,7000)
        ax1.set_ylim(-0.04,0.04)
        ax1.set_xlabel("Time(sec)",fontsize=15)
        ax1.set_ylabel("Signal (m)",fontsize=15)
        ax2.set_xlabel("frequency(Hz)",fontsize=15)
        ax2.set_ylabel("Power(normed)",fontsize=15)
        ax2.legend(loc='best')
        ax1.set_title('Signal vs. time, excited at %.2f'%self.x0,fontsize=18)
        ax2.set_title('power spectrum',fontsize=18)        
        plt.show() 
    def cal(self):
        self.y1=[0]*(self.M+1)
        for i in range(1,self.M,1):
            self.y1[i]=0.05*exp(-(i*self.dx-0.4)**2/0.05**2)
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
        self.ax.set_title('Wave on a string: Triangle',fontsize=18)
        self.ax.legend(loc='best',fontsize=14)
        plt.savefig(r'd:\save'+'%.2f'%_t+r'.jpg')

cmp=WAVE(0.5)
cmp.cal_spec(0.5)
cmp.plot_spec()     
        