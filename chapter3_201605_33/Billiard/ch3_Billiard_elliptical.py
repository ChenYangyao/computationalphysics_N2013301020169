from numpy import * 
import matplotlib.pyplot as plt

# for a elliptical-bounded boundary
class BILLIARD(object):
    def __init__(self,_x0=sqrt(2),_y0=0.,_vx0=0,_vy0=1,_dt=0.001,_time=500):
        self.dt, self.time, self.n = _dt, _time, int(_time/_dt)
        self.x, self.y, self.vx, self.vy = [_x0], [_y0], [_vx0], [_vy0]
    def cal(self):
        for i in range(self.n):
            self.nextx = self.x[-1]+self.vx[-1]*self.dt
            self.nexty = self.y[-1]+self.vy[-1]*self.dt
            self.nextvx, self.nextvy = self.vx[-1], self.vy[-1]
            if self.nextx**2/3+self.nexty**2>1:
                self.nextx=self.x[-1]
                self.nexty=self.y[-1]
                while not(self.nextx**2/3+self.nexty**2>1):
                    self.nextx=self.nextx+self.nextvx*self.dt/100
                    self.nexty=self.nexty+self.nextvy*self.dt/100
                self.norm=array([self.nextx,3*self.nexty])
                self.norm=self.norm/sqrt(dot(self.norm,self.norm))
                self.v=array([self.nextvx,self.nextvy])
                self.v_perpendicular=dot(self.v,self.norm)*self.norm
                self.v_parrallel=self.v-self.v_perpendicular
                self.v_perpendicular=-self.v_perpendicular
                [self.nextvx,self.nextvy]=self.v_parrallel+self.v_perpendicular
            self.x.append(self.nextx)
            self.y.append(self.nexty)
            self.vx.append(self.nextvx)
            self.vy.append(self.nextvy)
    def plot_position(self,_ax):
        _ax.plot(self.x,self.y,'-b',label='Ellipse'+r'$a=2,b=1$')
        _ax.plot(sqrt(3)*cos(linspace(0,2*pi,200)),sin(linspace(0,2*pi,200)),'-k',lw=5)
    def plot_phase(self,_ax,_secy=0):
        self.secy=_secy
        self.phase_x, self.phase_vx = [], []
        for i in range(len(self.x)):
            if abs(self.y[i]-self.secy)<1E-3 and abs(self.vx[i])<0.95:
                self.phase_x.append(self.x[i])
                self.phase_vx.append(self.vx[i])
        _ax.plot(self.phase_x,self.phase_vx,'ob',markersize=2,label='Ellipse'+r'$a=2,b=1$')
        
fig= plt.figure(figsize=(10,4))
ax1=plt.axes([0.1,0.15,0.35,0.7])
ax2=plt.axes([0.6,0.15,0.35,0.7])
ax1.set_xlim(-2,2)
ax1.set_ylim(-1.1,1.1)
ax2.set_xlim(-2,2)
ax2.set_ylim(-1.1,1.1)  
ax1.set_title('Elliptical boundary: '+r'$a^2 =3,b^2 = 1$',fontsize=18)
ax2.set_title('Phase-space: '+r'$a^2 =3,b^2 = 1$',fontsize=18)
ax1.set_xlabel(r'$x(m)$',fontsize=18)
ax1.set_ylabel(r'$y(m)$',fontsize=18)
ax2.set_xlabel(r'$x(m)$',fontsize=18)
ax2.set_ylabel(r'$v_x (m/s)$',fontsize=18)
cmp=BILLIARD()
cmp.cal()
cmp.plot_position(ax1)
cmp.plot_phase(ax2)
plt.show()