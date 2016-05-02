'''
PROGRAM rigidball
Author: Chen Yangyao     Last Modify:20160502
this program solves the rigid model for molecules
''' 
from numpy import * 
import matplotlib.pyplot as plt

# class RIGID : solves for tje rigid model of molecules
# where:  
#         number : number of molecule in a line 
#         radii : radii of molecule 
#         v0,sigma : initial velocity and its distribution variance
class RIGID(object):
    def __init__(self,_number=15,_radii=0.02,_v0=0.3,_sigma=0.1):
        self.n=_number**2
        self.r=_radii
        self.v=array([])
        self.theta=array([])
        for i in range(self.n):
            self.v=hstack((self.v,random.normal(_v0,_sigma)))
            self.theta=hstack((self.theta,2*pi*random.random()))
        self.vx=[self.v*cos(self.theta)]
        self.vy=[self.v*sin(self.theta)]
        self.x=[array([])]
        self.y=[array([])]
        for i in linspace(-0.3,0.3,_number):
            for j in linspace(-0.3,0.3,_number):
                self.x[0]=hstack((self.x[0],i))
                self.y[0]=hstack((self.y[0],j))
    def cal(self,_dt=0.1,_time=100):         # Use Euler method to solve the rigid model 
        self.dt=_dt
        self.time=_time
        for k in range(int(self.time/self.dt)):
            print k,'in',int(self.time/self.dt)
            self.nextx=self.x[-1]+self.vx[-1]*self.dt
            self.nexty=self.y[-1]+self.vy[-1]*self.dt
            self.nextvx=self.vx[-1]
            self.nextvy=self.vy[-1]
            for i in range(self.n):
                for j in range(i+1,self.n,1):
                    if ((self.nextx[i]-self.nextx[j])**2+(self.nexty[i]-self.nexty[j])**2 < (2*self.r)**2):
                        self.norm = array([self.nextx[i]-self.nextx[j],self.nexty[i]-self.nexty[j]])
                        self.norm=self.norm/sqrt(dot(self.norm,self.norm))
                        self.vi=array([self.nextvx[i],self.nextvy[i]])
                        self.vj=array([self.nextvx[j],self.nextvy[j]])
                        self.vi_perp = dot(self.vi,self.norm)*self.norm
                        self.vi_para = self.vi- self.vi_perp
                        self.vj_perp = dot(self.vj,self.norm)*self.norm
                        self.vj_para = self.vj- self.vj_perp
                        self.vi=self.vi_para-self.vi_perp
                        self.vj=self.vj_para-self.vj_perp
                        [self.nextvx[i],self.nextvy[i]]=self.vi
                        [self.nextvx[j],self.nextvy[j]]=self.vj
                        self.nextx[i]=self.x[-1][i]
                        self.nexty[i]=self.y[-1][i]
                        self.nextx[j]=self.x[-1][j]
                        self.nexty[j]=self.y[-1][j]
            for i in range(self.n):
                if abs(self.nextx[i])>(0.5-self.r):
                    self.nextx[i]=self.nextx[i]-sign(self.nextx[i])*(1-self.r*2)
                if abs(self.nexty[i])>(0.5-self.r):
                    self.nexty[i]=self.nexty[i]-sign(self.nexty[i])*(1-self.r*2)
            self.gama=sqrt(self.nextvx**2+self.nextvy**2)
            self.gama=abs(self.gama-(3*1E-4))/self.gama        # to lower the temperature
            self.nextvx=self.nextvx*self.gama
            self.nextvy=self.nextvy*self.gama
            self.x.append(self.nextx)
            self.y.append(self.nexty)
            self.vx.append(self.nextvx)
            self.vy.append(self.nextvy)
    def plot(self,_ax,_i):         # plot the ith picture
        _ax.plot(self.x[_i],self.y[_i],'ob',markersize=4,label='last second')

cmp=RIGID()
cmp.cal()
'''
for i in range(5):
    fig=plt.figure(figsize=(6,6))
    ax1=plt.subplot(111)        
    cmp.plot(ax1,i*1000)
    ax1.set_xlim(-0.7,0.7)
    ax1.set_ylim(-0.7,0.7)
    plt.savefig("C:\\Users\\ChenYangyao\\Desktop\\pic_ch3_lorenz\\rigid\\i_%d .png"%i)
'''
# compute and output data
myfile=file(r'd:\datax.txt','w')
for i in range(len(cmp.x)):
    for j in range(cmp.n):
        print >> myfile,cmp.x[i][j]
myfile.close()
myfile=file(r'd:\datay.txt','w')
for i in range(len(cmp.y)):
    for j in range(cmp.n):
        print >> myfile,cmp.y[i][j]
myfile.close()
                    
                    
        
        