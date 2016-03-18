# PROGRAM radioactive_decay_of_two_type_nuclei

# Purpose: This program solves for the decay of two type nuclei.
#     The number of two type nuclei is detremined by two coupled ODEs (with both first order):
#                     dn_A/dt = (n_B - n_A)/tau 
#                     dn_B/dt = (n_A - n_B)/tau
#     where tau is decay constant
#     We solve the ODEs with Euler method. 
     
# Author: Chen Yangyao    
# Last modify: 20160318

# import the tool
import pickle
from pylab import *

# declare the var. & para.
N_A = []    
N_B = []    # N_A & N_B are number of two type nuclei
t = []      # time
tau = 0
time = 0
n = 1
dt= 1        # tau: decay constant, dt: numerical step size, time: total time, n:total number of steps   
N_A_total = []   # storage of all data 
N_B_total = []
t_total = []
dt_total = []
N_A_diff, N_B_diff = [],[]
            # store the difference between theory & simulation
times = int(raw_input('times of differnt dt (no less than 1) :  '))
            # dt ranges
color = ['ob','xg','*y','oy','pk','+m','xc']+['sy']*(times-7)
            # color sets of plot
# to set zero of some var. & para. 
def delete(): 
    global dt,n,t,N_A,N_B
    dt = 1.
    n = 1.
    N_A = [N_A[0]]
    N_B = [N_B[0]]
    t = [t[0]]
    
    
# initialize these var. & para.
def initialize(N_A,N_B,t,_tau,_dt,_time,_n): 
    global tau,dt,time,n
    N_A.append(float(raw_input('initial number of nuclei(type A):  ')))
    N_B.append(float(raw_input('initial number of nuclei(type B):  ')))
    t.append(float(raw_input('initial time:    ')))
    tau = float(raw_input('decay constant:    '))
    time = float(raw_input('numerical total time:  '))

# calculate the time dependence of N_A & N_B using Euler method
def calculate(N_A,N_B,t,tau,dt,time,n):
    for i in range(n): 
        N_A.append(N_A[-1]+(N_B[-1]-N_A[-1])/tau*dt)
        N_B.append(N_B[-1]+(N_A[-1]-N_B[-1])/tau*dt)
        t.append(t[-1]+dt)

# store data:
#      pickle data, stored in 'nuclei_decay_two_type_pickle.txt'
#      txt data, stored in 'nuclei_decay_two_type_txt.txt'
def store(N_A,N_B,t,tau,dt,time,n):
    myfile = open('nuclei_decay_two_type_pickle.txt','a')
    pickle.dump(t,myfile)
    pickle.dump(N_A,myfile)
    pickle.dump(N_B,myfile)
    myfile.close()
    myfile = open('nuclei_decay_two_type_txt.txt','a')
    print >> myfile,'This data is a numerical solution for the decay problem'
    for i in range(len(t)): 
        print >> myfile,t[i],N_A[i],N_B[i]


# MAIN: call the function & do all loops
# compute number of radioactive nuclei
initialize(N_A,N_B,t,tau,dt,time,n)
c1 = (N_B[0]+N_A[0])/2.
c2 = (N_A[0]-N_B[0])/2.
for i in range(times): 
    delete()
    dt = float(raw_input('numerical step size:  '))
    n = int((time/dt))
    calculate(N_A,N_B,t,tau,dt,time,n)
    store(N_A,N_B,t,tau,dt,time,n)
    N_A_total.append(N_A)
    N_B_total.append(N_B)
    t_total.append(t)
    dt_total.append(dt)
    N_A_theory = c2*np.exp(-2/tau*(array(t)-t[0]))+c1
    N_B_theory = -c2*np.exp(-2/tau*(array(t)-t[0]))+c1
    N_A_diff.append(array(N_A_total[i])-array(N_A_theory))
    N_B_diff.append(array(N_B_total[i])-array(N_B_theory))

# plot graphics to show the numerical result
fig_all=figure(figsize=(14,7),dpi=80)
# subplot1: plot time dependence of number of nuclei
subplot(121)
for i in range(times): 
    plot(t_total[i],N_A_total[i],color[i],label=r'$dt= $'+str(dt_total[i]))
    plot(t_total[i],N_B_total[i],color[i])
plot(t,N_A_theory,'-r',linewidth=2,label='theory')
plot(t,N_B_theory,'-r',linewidth=2)
xlim(t[0]-time/10.,t[0]+time+time/10)
ymin = min([N_A[0],N_B[0]])
ymax = max([N_A[0],N_B[0]])
ylim(ymin-(ymax-ymin)/10,ymax+(ymax-ymin)/10)
xlabel('time')
ylabel('number of nuclei')
title('decay: two types of nucleus',fontsize=20)
legend(loc = 'best')
annotate(r'$N_A$',xy=(t[int(len(t)/2)],
         N_A_theory[int(len(N_A_theory)/2)]), 
         xytext=(+30, +100),
         textcoords='offset points', fontsize=16,
         arrowprops=dict(facecolor='black',shrink=0.1))
annotate(r'$N_B$',xy=(t[int(len(t)/2)],
         N_B_theory[int(len(N_B_theory)/2)]), 
         xytext=(+30, -100),
         textcoords='offset points', fontsize=16,
         arrowprops=dict(facecolor='black',shrink=0.1))
subplot(122)
 # subplot2: show the difference between theory and simulation
for i in range(times): 
    plot(t_total[i],N_A_diff[i],color[i],label=r'$dt= $'+str(dt_total[i]))
    plot(t_total[i],N_B_diff[i],color[i])
xlim(t[0]-time/10.,t[0]+time+time/10)
ymin=min([min(N_A_diff[0]),min(N_A_diff[-1]),min(N_B_diff[0]),min(N_B_diff[-1])])
ymax=max([max(N_A_diff[0]),max(N_A_diff[-1]),max(N_B_diff[0]),max(N_B_diff[-1])])
ylim(ymin-(ymax-ymin)/10,ymax+(ymax-ymin)/10)
xlabel('time')
ylabel('difference')
title('difference between theory and simulation',fontsize=20)
legend(loc ='best')

savefig('decay of nuclei.jpg')
show(fig_all)
