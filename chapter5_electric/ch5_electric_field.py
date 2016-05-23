'''
program : electric field

'''
from numpy import *
import mpl_toolkits.mplot3d
import matplotlib.pyplot as plt

class ELECTRIC_FIELD(object):
    def __init__(self, V1, V2, n):
        self.V1=float(V1)
        self.V2=float(V2)
        self.n=n
        self.V=[]
        for j in range(self.n):
            self.V.append([self.V1]+[0.]*(self.n-2)+[self.V2])
        for i in range(self.n):
            self.V[0][i]=self.V1+(self.V2-self.V1)/(self.n-1)*i
            self.V[self.n-1][i]=self.V1+(self.V2-self.V1)/(n-1)*i
    def update_V(self):
        self.alpha=2./(1.+pi/self.n)
        self.counter=0
        while True:
            self.delta_V=0.
            for j in range(1,self.n-1):
                for i in range(1,self.n-1):
                    self.next_V=1./4.*(self.V[j-1][i]+self.V[j+1][i]+self.V[j][i-1]+self.V[j][i+1])
                    self.V[j][i]=self.alpha*(self.next_V-self.V[j][i])+self.V[j][i]
                    self.delta_V=self.delta_V+abs(self.alpha*(self.next_V-self.V[j][i]))
            self.counter=self.counter+1
            if (self.delta_V < abs(self.V2-self.V1)*(1.0E-5)*self.n*self.n and self.counter >= 10):
                break
    def plot_3d(self,ax,x1,x2,y1,y2):
        for j in range (self.n):
            for i in range(self.n):
                self.x=x1+(x2-x1)/(self.n-1.)*i
                self.y=y2-(y2-y1)/(self.n-1.)*j
                ax.plot([self.x],[self.y],[self.V[j][i]],'ok',markersize=2)
    def export_data(self,x1,x2,y1,y2):
        self.mfile=open(r'd:\data.txt','w')
        for j in range(self.n):
            for i in range(self.n):
                self.x=x1+(x2-x1)/(self.n-1.)*i
                self.y=y2-(y2-y1)/(self.n-1.)*j
                print >> self.mfile, self.x,self.y,self.V[j][i]
        self.mfile.close()
            
        
           
cmp=ELECTRIC_FIELD(-1.,1.,100)
cmp.update_V()
cmp.export_data(-1.,1.,-1.,1.)

fig=plt.figure(figsize=(7,7))
ax1= plt.subplot(1,1,1,projection='3d')
cmp.plot_3d(ax1,-1.,1.,-1.,1.)
plt.show(fig)

     
        