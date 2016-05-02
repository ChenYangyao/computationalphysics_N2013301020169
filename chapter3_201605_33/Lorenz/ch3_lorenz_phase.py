'''
PROGRAM lorenz_phase
Author: Chen Yangyao     Last Modify:20160502
this program solves the Lorenz model and give trajectory and phase-space plot
''' 

from numpy import * 
import matplotlib.pyplot as plt

#class:LORENZ solves for Lorenz model 
#where:   
#       r : driving    dt: time step 
#       time :total time 
#       sigma, b : other parameter    
class LORENZ(object):
    def __init__(self,_r=25,_dt=0.001,_time=50,_sigma=10,_b=10/3):
        self.t = [0]
        self.x, self.y, self.z = [1], [0], [0]
        self.r, self.sigma, self.b = _r, _sigma, _b
        self.dt, self.time, self.n = _dt, _time, int(_time/_dt)
        print 'Lorenz init, dt=%.3f, time=%.3f, r=%.3f '%(self.dt,self.time,self.r)
    def cal(self):           # use Euler method solves Lorenz Equations numerically
        for i in range(self.n):
            self.t.append(self.t[-1]+self.dt)           
            self.x.append(self.x[-1]+self.sigma*(self.y[-1]-self.x[-1])*self.dt)
            self.y.append(self.y[-1]+(-self.x[-2]*self.z[-1]+self.r*self.x[-2]-self.y[-1])*self.dt)
            self.z.append(self.z[-1]+(self.x[-2]*self.y[-2]-self.b*self.z[-1])*self.dt)
        print 'Lorenz cal',1
    def plot_z(self,_ax):           # plot trajectory z vs time
        _ax.plot(self.t,self.z,'-r',label='r= %.2f'%self.r)
        print 'z_t plot',1
    def plot_x(self,_ax):           # plot trajectory x vs time
        _ax.plot(self.t,self.x,'-r',label='r= %.2f'%self.r)
        print 'x_t plot',1
    def plot_y(self,_ax):           # plot trajectory y vs time
        _ax.plot(self.t,self.y,'-r',label='r= %.2f'%self.r)
        print 'y_t plot',1
    def plot_phase_zx(self,_ax):           # plot phase-space z vs x
        _ax.plot(self.x,self.z,'ob',markersize=0.3,label='r= %.2f'%self.r)
        print 'zx phase plot',1
    def plot_phase_xy(self,_ax):           # plot phase-space x vs y
        _ax.plot(self.y,self.x,'ob',markersize=0.3,label='r= %.2f'%self.r)
        print 'xy phase plot',1
    def plot_phase_yz(self,_ax):           # plot phase-space y vs z
        _ax.plot(self.z,self.y,'ob',markersize=0.3,label='r= %.2f'%self.r)
        print 'yz phase plot',1
    def plot_section_zx(self,_ax,_secy=0.):           # plot phase-space z vs x at y=0 section
        self.z_sec, self.x_sec = [], []
        self.secy = _secy
        for i in range(len(self.t)):
            if abs(self.y[i]-self.secy)<4E-3:
                self.z_sec.append(self.z[i])
                self.x_sec.append(self.x[i])
        if len(self.z_sec)>0:
            _ax.plot(self.x_sec,self.z_sec,'ok',markersize=4,label='r= %.2f'%self.r)
            print 'section plot,secy=',self.secy
            return 1
        else:
            print 'sec_y overranging'
            return 0
'''
# set r = 5.,15.,29.
# plot z versus time
fig= plt.figure(figsize=(8,8))
ax1=plt.subplot(311)
ax2=plt.subplot(312)
ax3=plt.subplot(313)
cmp=LORENZ(5.)
cmp.cal()
cmp.plot_z(ax1)
cmp=LORENZ(15.)
cmp.cal()
cmp.plot_z(ax2)
cmp=LORENZ(29.)
cmp.cal()
cmp.plot_z(ax3)
ax1.legend(loc='upper right')
ax2.legend(loc='upper right')
ax3.legend(loc='upper right')
ax1.set_title('Lorenz Model: z versus time',fontsize=18)
ax3.set_xlabel('time ',fontsize=18)
ax2.set_ylabel('z ',fontsize=18)
plt.show(fig)
'''

# plot phase space and section
fig= plt.figure(figsize=(10,10))
ax1=plt.axes([0.1,0.55,0.35,0.35])
ax2=plt.axes([0.6,0.55,0.35,0.35])
ax3=plt.axes([0.1,0.1,0.35,0.35])
ax4=plt.axes([0.6,0.1,0.35,0.35])

cmp=LORENZ(29.)
cmp.cal()
cmp.plot_phase_xy(ax1)
cmp.plot_phase_yz(ax2)
cmp.plot_phase_zx(ax3)
cmp=LORENZ(29.,0.001,1000)
cmp.cal()
cmp.plot_section_zx(ax4)

ax1.legend(loc='lower right')
ax2.legend(loc='lower right')
ax3.legend(loc='lower right')
ax4.legend(loc='lower right')
ax1.set_xlabel('y',fontsize=18)
ax2.set_xlabel('z',fontsize=18)
ax3.set_xlabel('x',fontsize=18)
ax4.set_xlabel('x',fontsize=18)
ax1.set_ylabel('x',fontsize=18)
ax2.set_ylabel('y',fontsize=18)
ax3.set_ylabel('z',fontsize=18)
ax4.set_ylabel('z',fontsize=18)
ax1.set_title('Phase space :x versus y',fontsize=18)
ax2.set_title('Phase space :y versus z',fontsize=18)
ax3.set_title('Phase space :z versus x',fontsize=18)
ax4.set_title('Section y=0 :z versus x',fontsize=18)

plt.show(fig)


            
'''
#  change r(driving force)
#  plot z vs t
for r in range(30):
    fig= plt.figure(figsize=(10,6))
    ax1=plt.subplot(111)
    cmp=LORENZ(r+1.)
    cmp.cal()
    cmp.plot_z(ax1)
    ax1.set_xlabel('t')
    ax1.set_ylabel('z')
    ax1.set_title('z versus t')
    plt.legend(loc='best')
    plt.savefig('C:\\Users\\ChenYangyao\\Desktop\\pic_ch3_lorenz\\z_t\\r_%.2f.png'%cmp.r)
'''



'''
#  change r(driving force)
#  plot phase zx, xy, yz
for r in range(20,30,1):
    fig= plt.figure(figsize=(10,6))
    ax1=plt.axes([0.1,0.6,0.3,0.3])
    ax2=plt.axes([0.6,0.1,0.3,0.3])
    ax3=plt.axes([0.1,0.1,0.3,0.3])
    cmp=LORENZ(r+1.,0.001,50)
    cmp.cal()
    cmp.plot_phase_zx(ax1)
    ax1.set_xlabel('x')
    ax1.set_ylabel('z')
    ax1.set_title('zx')
    cmp.plot_phase_xy(ax2)
    ax2.set_xlabel('y')
    ax2.set_ylabel('x')
    ax2.set_title('xy')
    cmp.plot_phase_yz(ax3)
    ax3.set_xlabel('z')
    ax3.set_ylabel('y')
    ax3.set_title('yz')
    ax1.legend(loc='best')
    ax2.legend(loc='best')
    ax3.legend(loc='best')
    plt.savefig('C:\\Users\\ChenYangyao\\Desktop\\pic_ch3_lorenz\\phase\\r_%.2f.png'%cmp.r)
'''



'''
#  change r(driving force)
#  plot section z vs x
for r in linspace(25,50,100):
    fig= plt.figure(figsize=(10,6))
    ax1=plt.subplot(111)
    cmp=LORENZ(r,0.001,300)
    cmp.cal()
    cmp.plot_section_zx(ax1,0.)
    ax1.set_xlabel('x')
    ax1.set_ylabel('z')
    ax1.set_title('zx')
    ax1.set_xlim(-30,30)
    ax1.set_ylim(0,80)
    ax1.legend(loc='best')
    plt.savefig('C:\\Users\\ChenYangyao\\Desktop\\pic_ch3_lorenz\\section\\r_%.2f.png'%r)
'''



'''
#  change y(driving force)
#  plot section z vs x
cmp=LORENZ(28,0.001,1000)
cmp.cal()
for secy in linspace(-30,30,60):
    fig= plt.figure(figsize=(10,6))
    ax1=plt.subplot(111)
    cmp.plot_section_zx(ax1,secy)
    ax1.set_xlabel('x')
    ax1.set_ylabel('z')
    ax1.set_title('zx')
    ax1.set_xlim(-15,15)
    ax1.set_ylim(5,40)
    ax1.legend(loc='best')
    plt.savefig('C:\\Users\\ChenYangyao\\Desktop\\pic_ch3_lorenz\\butter\\secy_%.2f.png'%secy)
'''
 