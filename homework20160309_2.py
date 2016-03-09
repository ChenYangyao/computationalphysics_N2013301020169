# program: it will give a plot of a sine curve
# the form of this curve is 15*sin(2*pi/30*x), it will move toward -x axis
# author: Chen Yangyao(2013301020169)    update time: 20160309

#----------- import the package --------
from math import *
import time
import os

#----------- set initial conditions ----
main_tab = [' '*75]*40    # the matrix that will be print
ii,jj,kk= 1,1,1           # loop var.

#----------- main program --------------
while kk<= 400:
    while jj<=75:         # this loop make the matrix like a sine curve
        fx = int(15*sin(2*3.14/30*jj+2*3.14/100*kk))
        tab_temp= main_tab[20-fx-1]
        tab_temp= list(tab_temp)
        tab_temp[jj-1]= '#'
        tab_temp= ''.join(tab_temp)
        main_tab[20-fx-1] = tab_temp
        jj= jj+1
    while ii<=40:         # print the matrix
        print main_tab[ii-1]
        ii=ii+1
    time.sleep(0.01)      
    os.system('cls')      # clear the screen
    ii,jj,kk=1,1,kk+1     # reset the initial conditions
    main_tab = [' '*75]*40
