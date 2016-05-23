'''
program : electric field near capacitor

'''

from numpy import *
import mpl_toolkits.mplot3d
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import time

class ELECTRIC_FIELD(object):
    def __init__(self, V1=1, V2=-1, V_boundary=0, n=30):
        self.V1=float(V1)
        self.V2=float(V2)
        self.V_boundary=float(V_boundary)
        self.n=int(n)
        self.s1, self.s3=int(n/3), int(n/3)
        self.s2, self.s4=int(self.n-2-2*self.s1), int(self.n-2*self.s3)
        self.V=[]
        for j in range(self.n):
            self.V.append([0.0]*self.n)
        for j in [0,self.n-1]:
            self.V[j]=[self.V_boundary]*self.n
        for j in range(self.s3,self.s3+self.s4):
            self.V[j][self.s1]=self.V1
            self.V[j][self.s1+self.s2+1]=self.V2            
    def update_V_Jacobi(self):
        self.counter=0        
        while True:
            self.V_next=[]
            for j in range(self.n):
                self.V_next.append([0.0]*self.n)
            for j in range(self.n):
                for i in range(self.n):
                    self.V_next[j][i]=self.V[j][i]
            self.delta_V=0.
            for j in range(1,self.n-1):
                for i in range(1,self.n-1):
                    if (j in range(self.s3,self.s3+self.s4)) and (i in [self.s1,self.s1+self.s2+1]):
                        continue
                    self.V_next[j][i]=1./4.*(self.V[j-1][i]+self.V[j+1][i]+self.V[j][i-1]+self.V[j][i+1])
                    self.delta_V=self.delta_V+abs(self.V_next[j][i]-self.V[j][i])
            self.counter=self.counter+1
            for j in range(self.n):
                for i in range(self.n):
                    self.V[j][i]=self.V_next[j][i]
            if (self.delta_V < abs(self.V2-self.V1)*(1.0E-5)*self.n*self.n and self.counter >= 10):
                break
        print 'jacobi itertion length n=',self.n,'  ',self.counter,' times'
        return self.counter
    def update_V_Gauss(self):
        self.counter=0
        while True:
            self.delta_V=0.
            for j in range(1,self.n-1):
                for i in range(1,self.n-1):
                    if (j in range(self.s3,self.s3+self.s4)) and (i in [self.s1,self.s1+self.s2+1]):
                        continue
                    self.next_V=1./4.*(self.V[j-1][i]+self.V[j+1][i]+self.V[j][i-1]+self.V[j][i+1])
                    self.delta_V=self.delta_V+abs(self.next_V-self.V[j][i])
                    self.V[j][i]=self.next_V   
            self.counter=self.counter+1
            if (self.delta_V < abs(self.V2-self.V1)*(1.0E-5)*self.n*self.n and self.counter >= 10):
                break
        print 'gauss itertion length n=',self.n,'  ',self.counter,' times'
        return self.counter    
    def update_V_SOR(self):
        self.alpha=2./(1.+pi/self.n)
        self.counter=0
        while True:
            self.delta_V=0.
            for j in range(1,self.n-1):
                for i in range(1,self.n-1):
                    if (j in range(self.s3,self.s3+self.s4)) and (i in [self.s1,self.s1+self.s2+1]):
                        continue
                    self.next_V=1./4.*(self.V[j-1][i]+self.V[j+1][i]+self.V[j][i-1]+self.V[j][i+1])
                    self.V[j][i]=self.alpha*(self.next_V-self.V[j][i])+self.V[j][i]
                    self.delta_V=self.delta_V+abs(self.alpha*(self.next_V-self.V[j][i]))
            self.counter=self.counter+1
            if (self.delta_V < abs(self.V2-self.V1)*(1.0E-5)*self.n*self.n and self.counter >= 10):
                break
        print 'SOR itertion length n=',self.n,'  ',self.counter,' times'
        return self.counter
    def Ele_field(self,x1,x2,y1,y2):
        self.dx=abs(x1-x2)/float(self.n-1)
        self.Ex=[]
        self.Ey=[]
        for j in range(self.n):
            self.Ex.append([0.0]*self.n)
            self.Ey.append([0.0]*self.n)
        for j in range(1,self.n-1,1):
            for i in range(1,self.n-1,1):
                self.Ex[j][i]=-(self.V[j][i+1]-self.V[j][i-1])/(2*self.dx)
                self.Ey[j][i]=-(self.V[j-1][i]-self.V[j+1][i])/(2*self.dx)
    def plot_3d(self,ax,x1,x2,y1,y2):
        self.x=linspace(x1,x2,self.n)
        self.y=linspace(y2,y1,self.n)
        self.x,self.y=meshgrid(self.x,self.y)
        self.surf=ax.plot_surface(self.x,self.y,self.V, rstride=1, cstride=1, cmap=cm.coolwarm)
        ax.set_xlim(x1,x2)
        ax.set_ylim(y1,y2)
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.01f'))
        ax.set_xlabel('x (m)',fontsize=14)
        ax.set_ylabel('y (m)',fontsize=14)
        ax.set_zlabel('Electric potential (V)',fontsize=14)
        ax.set_title('Potential near capacitor',fontsize=18)
    def plot_2d(self,ax1,ax2,x1,x2,y1,y2):
        self.x=linspace(x1,x2,self.n)
        self.y=linspace(y2,y1,self.n)
        self.x,self.y=meshgrid(self.x,self.y)
        
        cs=ax1.contour(self.x,self.y,self.V)
        plt.clabel(cs, inline=1, fontsize=10)
        ax1.set_title('Equipotential lines',fontsize=18)
        ax1.set_xlabel('x (m)',fontsize=14)
        ax1.set_ylabel('y (m)',fontsize=14)
        
        for j in range(1,self.n-1,1):
            for i in range(1,self.n-1,1):
                ax2.arrow(self.x[j][i],self.y[j][i],self.Ex[j][i]/40,self.Ey[j][i]/40,fc='k', ec='k')             
        ax2.set_xlim(-1.,1.)
        ax2.set_ylim(-1.,1.)
        ax2.set_title('Electric field',fontsize=18)
        ax2.set_xlabel('x (m)',fontsize=14)
        ax2.set_ylabel('y (m)',fontsize=14)
    def export_data(self,x1,x2,y1,y2):
        self.mfile=open(r'd:\data.txt','w')
        self.x=linspace(x1,x2,self.n)
        self.y=linspace(y2,y1,self.n)
        self.x,self.y=meshgrid(self.x,self.y)
        for j in range(self.n):
            for i in range(self.n):
                print >> self.mfile, self.x[j][i],self.y[j][i],self.V[j][i]
        self.mfile.close()
        
    
n_min=10
n_max=50
n_jacobi=[]
iters_jacobi=[]
time_jacobi=[]
for i in range(n_min,n_max,2):
    start=time.clock()
    cmp=ELECTRIC_FIELD(1,-1,0,i)
    iters_jacobi.append(cmp.update_V_Jacobi())
    end=time.clock()
    
    time_jacobi.append(end-start)   
    n_jacobi.append(i)


n_gauss=[]
iters_gauss=[]
time_gauss=[]
for i in range(n_min,n_max,2):
    start=time.clock()
    cmp=ELECTRIC_FIELD(1,-1,0,i)
    iters_gauss.append(cmp.update_V_Gauss())
    end=time.clock()
   
    time_gauss.append(end-start)   
    n_gauss.append(i)
    
n_SOR=[]
iters_SOR=[]
time_SOR=[]
for i in range(n_min,n_max,2):
    start=time.clock()
    cmp=ELECTRIC_FIELD(1,-1,0,i)
    iters_SOR.append(cmp.update_V_SOR())
    end=time.clock()
    time_SOR.append(end-start)
    n_SOR.append(i)
print time_SOR

mfile=open(r'd:\data.txt','w')
for i in range(len(n_jacobi)):
    print >> mfile, n_jacobi[i], time_jacobi[i], time_gauss[i], time_SOR[i]
mfile.close()


fig=plt.figure(figsize=(14,7))
ax1=plt.axes([0.1,0.1,0.35,0.8])
ax2=plt.axes([0.55,0.1,0.35,0.8])

ax1.plot(n_jacobi,iters_jacobi,'or',markersize=5)
ax1.plot(n_gauss,iters_gauss,'ob',markersize=5)
ax1.plot(n_SOR,iters_SOR,'oy',markersize=5)
ax2.plot(n_jacobi,time_jacobi,'or',markersize=5)
ax2.plot(n_gauss,time_gauss,'ob',markersize=5)
ax2.plot(n_SOR,time_SOR,'oy',markersize=5)

plt.show(fig)






        
        