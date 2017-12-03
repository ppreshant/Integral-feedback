# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 23:39:52 2017

@author: Prashant
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Rate constants
k = 10; k2 = 1.4513; k3 = 2.3679;
g1 = 1.2337; g2 = 3.0155; g3 = 1.1114;
du = 0;
mu = 2; al = .081; 


# Solving the ODEs for the reaction given the rate consts
# and initial conditions
def OdeReact(t0,t1,X0, k,k2,k3,g1,g2,g3,du,mu,al):
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
    p = 0
    return time, MA, sp, p    
    
    
# first step

param = k,k2,k3,g1,g2,g3,du,mu,al
X0 = np.array([0, 0, 0, 0.1])
time1, MA1, sp1, p1 = OdeReact(0,100,X0,*param)


# second step

mu = 5
param = k,k2,k3,g1,g2,g3,du,mu,al
time2, MA2, sp2, p2 = OdeReact(100,150,MA1[-1],*param)


# Third step

mu = 1
param = k,k2,k3,g1,g2,g3,du,mu,al
time3, MA3, sp3, p3 = OdeReact(150,200,MA2[-1],*param)


# Pack all steps into 1 vector / matrix 

MAf = np.row_stack([MA1,MA2,MA3])
tf = np.append(np.append(time1,time2),time3)
spf = sp1+sp2+sp3

plt.plot(tf,MAf[:,-2:])
plt.ylabel('Y(t)')
plt.xlabel('Time (s)')
plt.plot(tf,spf,'--k', linewidth=1)
plt.legend(['x3','V','Mu'],loc='upper right')

#plt.plot(tf,MAf)
#plt.ylabel('Y(t)')
#plt.xlabel('Time (s)')
#plt.plot(tf,spf,'--k', linewidth=1)
#plt.legend(['x1','x2','x3','V','Mu'],loc='upper right')


