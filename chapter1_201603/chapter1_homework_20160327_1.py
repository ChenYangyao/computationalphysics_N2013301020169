# PROGRAM radioactive_decay_of_two_type_nuclei

# Purpose: This program solves for the decay of two type nuclei.
#     The number of two type nuclei is detremined by two coupled ODEs (with both first order):
#                     dn_A/dt = (n_B - n_A)/tau 
#                     dn_B/dt = (n_A - n_B)/tau
#     where tau is decay constant
#     We solve the ODEs with Euler method. 

# Author: Chen Yangyao    
# Last modify: 20160327
import pylab as pl   # import the tools
import pickle
import numpy as np
class DECAY(object):
    '''this program solve for the radioactive decay 
	   two types of nuclei will be considered
    '''
    def __init__(self,na=1000., nb=0., time=30.,tau=10., dt=0.5):
        self.na = float(na)   # na,nb: initial number of nuclei
        self.nb = float(nb)
        self.tau = float(tau)  # tau: time constant
        self.dt = float(dt)    # dt: step size
        self.time = float(time)# time :total time
        self.n = int(time/dt)  # n: step number
    def calculate(self):       # na_t,nb_t: # storage of all data 
        self.na_t = [self.na]
        self.nb_t = [self.nb]
        self.t_t = [0.]
        for ii in range(self.n):
            self.na_t.append((self.nb_t[-1]-self.na_t[-1])/self.tau*self.dt+self.na_t[-1])
            self.nb_t.append((self.na_t[-1]-self.nb_t[-1])/self.tau*self.dt+self.nb_t[-1])
            self.t_t.append(self.t_t[-1]+self.dt)
    def store(self,dt):     # output all data
                            # it will in dir d:\
        mfile = open(r'd:\decay_twotype_%s.txt'%dt,'w')
        pickle.dump(self.t_t,mfile)
        pickle.dump(self.na_t,mfile)
        pickle.dump(self.nb_t,mfile)
        mfile.close()
    def theory(self):     
        self.c1 = (self.na+self.nb)/2.
        self.c2 = (self.na-self.nb)/2.
        self.na_theory = [na]   # na_theory, nb_theory: theory results
        self.nb_theory = [nb]
        for ii in range(self.n):
            self.na_theory.append(self.c2*np.exp(-2./self.tau*(ii+1)*self.dt)+self.c1)
            self.nb_theory.append(-self.c2*np.exp(-2./self.tau*(ii+1)*self.dt)+self.c1)
    def graph(self,cmd,dt):     # plot number of nuclei
        pl.plot(self.t_t,self.na_t,cmd,label='d'+r'$t=%s$'%dt+'s')
        pl.plot(self.t_t,self.nb_t,cmd)
    def graph_theory(self,cmd):  # plot theory graphic
        pl.plot(self.t_t,self.na_theory,cmd,linewidth=2,label='Theory')
        pl.plot(self.t_t,self.nb_theory,cmd,linewidth=2)
    def graph_diff(self,cmd,dt):  # plot deviation of theory & discrete var.
        pl.plot(self.t_t,np.array(self.na_t)-np.array(self.na_theory),cmd,label='d'+r'$t=%s$'%dt+'s')
        pl.plot(self.t_t,np.array(self.nb_t)-np.array(self.nb_theory),cmd)
        
na = 1000.      # initialize all var. & para.
nb = 10.
time = 30.
tau = 10.
dt = [5.,2.,1.,2.,0.5]   # different step sizes
cmd = ['ob','xg','*y','oy','pk','+m','xc']  # plot command for different lines
figure = pl.figure(figsize=(14,7),dpi=80)
for jj in range(5) :     # plot numerical data and deviation
    cal = DECAY(na,nb,time,tau,dt[jj])    
    cal.calculate()
    cal.store(dt[jj])
    cal.theory()
    pl.subplot(121)
    cal.graph(cmd[jj],dt[jj])
    pl.subplot(122)
    cal.graph_diff(cmd[jj],dt[jj])
pl.subplot(121)  # theoretical line
cal.graph_theory('r')
pl.xlabel('Time(s)',fontsize=20)
pl.ylabel('Number of nuclei',fontsize=20)
pl.title('decay of two types of nucleus',fontsize=25)
pl.legend(loc = 'best')
pl.annotate(r'$N_A$',xy=(13,550),    # an arrow in plot
         xytext=(+30, +50),
         textcoords='offset points', fontsize=16,
         arrowprops=dict(facecolor='black',shrink=0.1))
pl.annotate(r'$N_B$',xy=(15,480), 
         xytext=(+30, -80),
         textcoords='offset points', fontsize=16,
         arrowprops=dict(facecolor='black',shrink=0.1))    
pl.subplot(122)
pl.xlabel('Time(s)',fontsize=20)    
pl.ylabel('Numerical error number',fontsize=20)
pl.title('Numerical errors',fontsize=25)
pl.legend(loc ='best')      
pl.savefig(r'd:\decay of nuclei.jpg') 
pl.show(figure)       
        
