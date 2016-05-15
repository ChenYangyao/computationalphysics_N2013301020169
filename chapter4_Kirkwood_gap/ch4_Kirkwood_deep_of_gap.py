'''
PROGRAM Kirkwood trajectory
this program solves for the motion of asteroids near the Kirkwood gap
both gravity from the Sun and Jupiter will be considered
'''

        
from numpy import *
import matplotlib.pyplot as plt


class KIRKWOOD(object):
    '''
    class KIRKWOOD solves for the asteroids' motion near the KIRKWOOD gap
    where :
                r0: initial radius of asteroid
                dt, time: time step size and total time
    use Euler-Cromer method, and ignore the perturbation of Jupiter by asteroids
    '''
    def __init__(self, _r0=3.276, _time=200., _dt=0.005):
        self.x_j, self.y_j, self.vx_j, self.vy_j=[5.2], [0], [0], [sqrt(4*pi*pi/5.2)]
        self.x_a, self.y_a, self.vx_a, self.vy_a=[_r0],[0],[0],[sqrt(4*pi*pi/_r0)]
        self.dt=_dt
        self.time=_time
        self.n=int(_time/_dt)
        self.Mj=318./332900.
        self.r0=_r0
        print 'KIRKWOOD init r0= %.5f'%_r0,1
    def cal(self):
        '''
        function cal (use euler method) solve the motion of asteroid and Jupiter
        '''
        for i in range(self.n):
            self.r_j=sqrt(self.x_j[-1]**2+self.y_j[-1]**2)
            self.r_a=sqrt(self.x_a[-1]**2+self.y_a[-1]**2)
            self.r_aj=sqrt((self.x_a[-1]-self.x_j[-1])**2+(self.y_a[-1]-self.y_j[-1])**2)
            self.vx_j.append(self.vx_j[-1]-self.dt*(4*pi*pi/self.r_j**3)*self.x_j[-1])
            self.vy_j.append(self.vy_j[-1]-self.dt*(4*pi*pi/self.r_j**3)*self.y_j[-1])
            self.vx_a.append(self.vx_a[-1]-self.dt*(4*pi*pi/self.r_a**3)*self.x_a[-1]-self.dt*(4*pi*pi*self.Mj/self.r_aj**3)*(self.x_a[-1]-self.x_j[-1]))
            self.vy_a.append(self.vy_a[-1]-self.dt*(4*pi*pi/self.r_a**3)*self.y_a[-1]-self.dt*(4*pi*pi*self.Mj/self.r_aj**3)*(self.y_a[-1]-self.y_j[-1]))
            self.x_j.append(self.x_j[-1]+self.vx_j[-1]*self.dt)
            self.y_j.append(self.y_j[-1]+self.vy_j[-1]*self.dt)
            self.x_a.append(self.x_a[-1]+self.vx_a[-1]*self.dt)
            self.y_a.append(self.y_a[-1]+self.vy_a[-1]*self.dt)
        print 'KIRKWOOD calculate',1
    def plot_asteroid(self,_ax):
        '''
        function plot_asteroid plots the trajectory of asteroid
        '''
        self.x_a_plot,self.y_a_plot=[],[]       
        for i in range(int(len(self.x_a)/10)):
            self.x_a_plot.append(self.x_a[i*10])
            self.y_a_plot.append(self.y_a[i*10])
        _ax.plot(self.x_a_plot,self.y_a_plot,'o',markersize=0.1,label=r'$r_0= $'+'%.2f'%self.r0)
    def plot_jupiter(self,_ax):
        '''
        function plot_asteroid plots the trajectory of jupiter
        '''
        self.x_j_plot,self.y_j_plot=[],[]
        for i in range(int(len(self.x_a)/10)):
            self.x_j_plot.append(self.x_j[i*10])
            self.y_j_plot.append(self.y_j[i*10])
        _ax.plot(self.x_j_plot,self.y_j_plot,'o',markersize=0.1)
    def find_deviation(self):
        '''
        function find_deviation gives how much the orbit ot asteroid deviate from circular orbit
        '''
        self.r_array=sqrt(array(self.x_a)**2+array(self.y_a)**2)
        return (max(self.r_array)-min(self.r_array))/self.r_array[0]



fig=plt.figure(figsize=(12,6))
ax1=plt.axes([0.1,0.1,0.8,0.8])

radius=[]
deviation=[]
for i in linspace(2.80,2.817,200):
    temp=KIRKWOOD(i)
    temp.cal()
    deviation.append(temp.find_deviation())
    radius.append(i)
    
mfile=open('d:\datach4r.txt','w')
for i in range(len(radius)):
    print >> mfile,radius[i]
mfile.close()
mfile=open('d:\datach4d.txt','w')
for i in range(len(radius)):
    print >> mfile,deviation[i]
mfile.close()

ax1.plot(radius,deviation,'-p',markersize=3,label='gap')

ax1.set_title('KIRKWOOD gap',fontsize=18)
ax1.set_xlabel('r  /AU',fontsize=18)
ax1.set_ylabel(r'$\Delta r / r_0$',fontsize=18)

ax1.legend(loc='best',fontsize=18)

plt.show(fig)
