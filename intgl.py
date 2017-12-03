# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 15:40:15 2017

@author: Prashant
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Rate constants
k = 1; k2 = 1; k3 = 1;
g1 = 1; g2 = 3; g3 = 1;
du = 1; bv = 1; gv = 1.5
mu = 2; al = .1; 


# Solving the ODEs for the reaction given the rate consts
# and initial conditions
def OdeReact(t0,t1,X0, k,k2,k3,g1,g2,g3,du,mu,al,bv,gv):
    def deriv(X, t):
    
        Ab = np.array([[-g1,    0,        0,    k],  
                   [ k2, -(g2 + k3),  0,    0],
                   [ 0,     k3,      -g3,   0],
                   [ 0,     0,   -al*X[3],  al*mu]])
        c = np.array([du,    0,        0,    0])                 
        return np.dot(Ab, X) + c

    time = np.linspace(t0, t1, t1 - t0 + 1)
    
    MA = odeint(deriv, X0, time)
    sp = [mu for i in time] # set point
    p = [g3 for i in time]  # Other parameters changed
    return time, MA, sp, p    
    
    
# first step

param = k,k2,k3,g1,g2,g3,du,mu,al,bv,gv
X0 = np.array([0, 0, 0, 0.1])
time1, MA1, sp1, p1 = OdeReact(0,100,X0,*param)


# second step

#mu = 4
g3 = g3/10.0
param = k,k2,k3,g1,g2,g3,du,mu,al,bv,gv
time2, MA2, sp2, p2 = OdeReact(100,150,MA1[-1],*param)


# Third step

#mu = 1
param = k,k2,k3,g1,g2,g3,du,mu,al,bv,gv
time3, MA3, sp3, p3 = OdeReact(150,200,MA2[-1],*param)


# Pack all steps into 1 vector / matrix 

MAf = np.row_stack([MA1,MA2,MA3])
tf = np.append(np.append(time1,time2),time3)
spf = sp1+sp2+sp3
pf = p1 + p2 + p3

#plt.plot(tf,MAf[:,-2:-1]) # plot x3
#plt.ylabel('Y(t)')
#plt.xlabel('Time (s)')
#plt.plot(tf,spf,'--k', linewidth=1)
##plt.plot(tf,pf,':r', linewidth=2)
#plt.legend(['x3','Mu'],loc='upper right')

plt.plot(tf,MAf[:,-2:]) # plot v and x3
plt.ylabel('Y(t)')
plt.xlabel('Time (s)')
plt.plot(tf,spf,'--k', linewidth=1)
plt.plot(tf,pf,':r', linewidth=2)
plt.legend(['x3','V','Mu','g3'],loc='upper right')

#plt.plot(tf,MAf)
#plt.ylabel('Y(t)')
#plt.xlabel('Time (s)')
#plt.plot(tf,spf,'--k', linewidth=1)
#plt.legend(['x1','x2','x3','V','Mu'],loc='upper right')


