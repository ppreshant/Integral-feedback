# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 01:26:14 2017


@author: Prashant
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Rate constants
k = 1; k2 = 1; k3 = 1;
g1 = 2; g2 = 3; g3 = 1;
du = 1; 
# control parameters
n = 10
mu = 3; te = 1; 

tmul = 100 ; # time to equilibriate the disturbances
par = 'k2' # parameter of interest
dp = {'k': k,'k2': k2,'k3': k3,'g1': g1,'g2': g2,'g3': g3,'du' : du, 'mu' : mu,'te':te,'n' : n} # dictionary with all parameters

# Solving the ODEs for the reaction given the rate consts
# and initial conditions
def OdeReact(t0,t1,X0, dp,par):
    [k3,k2,te,g3,g2,g1,k,du,mu,n] = dp.values() # order here is important
    def deriv(X, t):
    
        Ab = np.array([[-g1,    0,    0,    k,     0],
                   [ k2, -(g2 + k3),  0,    0,     0],
                   [ 0,     k3,      -g3,   0,     0],
                   [ 0,     0,        0,    0,    -n*X[3]],
                   [ 0,     0,       te, -n*X[4],  0]])
        
        
        c = np.array([du, 0, 0, te*mu, 0])                 
        return np.dot(Ab, X) + c

    time = np.linspace(t0, t1, t1 - t0 + 1)
    
    MA = odeint(deriv, X0, time)
    sp = [mu for i in time] # set point
    p = [dp[par] for i in time]  # Other parameters changed
    return time, MA, sp, p    
    
    
# first step

X0 = np.array([0.1, 1, .1, .1, .1])
time1, MA1, sp1, p1 = OdeReact(0,tmul,X0,dp,par)


# second step

dp[par] = dp[par]*10
time2, MA2, sp2, p2 = OdeReact(tmul,2*tmul,MA1[-1],dp,par)


# Third step

dp[par] = dp[par]*.01
time3, MA3, sp3, p3 = OdeReact(2*tmul,3*tmul,MA2[-1],dp,par)


# Pack all steps into 1 vector / matrix 

MAf = np.row_stack([MA1,MA2,MA3])
tf = np.append(np.append(time1,time2),time3)
spf = sp1+sp2+sp3
pf = p1 + p2 + p3

plt.plot(tf,MAf[:,-3:-2]) # plot x3
plt.ylabel('Y(t)')
plt.xlabel('Time (s)')
plt.plot(tf,spf,'--k', linewidth=1)
plt.plot(tf,pf,':r', linewidth=2)
plt.legend(['x3','Mu',par],loc='best')

#plt.plot(tf,MAf[:,-3:]) # plot v and x3
#plt.ylabel('Y(t)')
#plt.xlabel('Time (s)')
#plt.plot(tf,spf,'--k', linewidth=1)
#plt.plot(tf,pf,':r', linewidth=2)
#plt.legend(['x3','Z1','Z2','Mu',par],loc='upper right')

#plt.plot(tf,MAf) # Plot everything
#plt.ylabel('Y(t)')
#plt.xlabel('Time (s)')
#plt.plot(tf,spf,'--k', linewidth=1)
#plt.legend(['x1','x2','x3','Z1','Z2','Mu'],loc='upper right')



