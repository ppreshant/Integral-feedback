# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 01:26:14 2017


@author: Prashant
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Rate constants : System parameters
du = 1.0; k2 = 1.0; k3 = 1.0;
g1 = 2.0; g2 = 3.0; g3 = 1.0;

# control parameters
k = .1; n = 1.0; mu = 3.0; theta = 10.0; # general
k4 = 0; k5 = 0; g4 = 0; g5 = 0; # synth and degradation of control species
km = 1.0; h = 1.0;  # hill fn

tmul = 200 ; # time to equilibriate the disturbances
par = 'k2' # parameter of interest
parscale = 40.0 # how much should the parameter be increased

dp = {'k': k,'k2': k2,'k3': k3,'g1': g1,'g2': g2,'g3': g3,'du' : du, 'mu' : mu,'theta':theta,'n' : n, 'km' : km,'h' : h, 'k4' : k4, 'k5' : k5, 'g4' : g4, 'g5' : g5} # dictionary with all parameters

# differential equation sets
def anti(X, t, *param):
    
#    [k3,k2,te,g3,g2,g1,k,du,mu,n] = param # order here is important : print( str(a)).translate(None,"'")
    [mu, g5, g4, g3, g2, g1, h, k, km, n, k3, k2, k5, k4, theta, du] = param;
    Ab = np.array([[-g1,    0,    0,    k,     0],
               [ k2, -(g2 + k3),  0,    0,     0],
               [ 0,     k3,      -g3,   0,     0],
               [ 0,     0,        0,    0,    -n*X[3]],
               [ 0,     0,       theta, -n*X[4],  0]])
    
    
    c = np.array([du, 0, 0, theta*mu, 0])                 
    return np.dot(Ab, X) + c

def hill(X, t, *param):
    
    [mu, g5, g4, g3, g2, g1, h, k, km, n, k3, k2, k5, k4, theta, du] = param;
    theta = 2.5;
    Ab = np.array([[-g1,    0,    0,    0,  0],  
               [ k2, -(g2 + k3),  0,    0,  0],
               [ 0,     k3,      -g3,   0,  0],
               [ 0,     0,   -theta*X[3],  -k4, 0],
               [0,  0,  0,  0,  0]])         
             
    c = np.array([du + k*pow(X[3]/km,h)/(1+pow(X[3]/km,h)) , 0, 0 , theta*mu,   0])      # constant production terms            
    return np.dot(Ab, X) + c

def intgl(X, t, *param):
    
    [mu, g5, g4, g3, g2, g1, h, k, km, n, k3, k2, k5, k4, theta, du] = param;
    theta = .04
    Ab = np.array([[-g1,    0,        0,    k,  0],  
               [ k2, -(g2 + k3),  0,    0,  0],
               [ 0,     k3,      -g3,   0,  0],
               [ 0,     0,   -theta*X[3],  theta*mu,    0],
               [0,0,0,0,0]])
    c = np.array([du,    0,        0,    0, 0])                 
    return np.dot(Ab, X) + c
    
# Solving the ODEs for the reaction given the rate consts
# and initial conditions
def OdeReact(deriv,t0,t1,X0, dpar,par):
    
    time = np.linspace(t0, t1, t1 - t0 + 1)
    
    MA = odeint(deriv, X0, time, tuple(dpar.values()))
    sp = [mu for i in time] # set point
    p = [dp[par] for i in time]  # Other parameters changed
    return time, MA, sp, p    
    
 
def quadsteps():    
    # first step
    
    X0 = np.array([0.1, 1, .1, .1, .1])
    time1, MA1, sp1, p1 = OdeReact(anti,0,tmul,X0,dp,par)
    
    
    # second step
    
    dp[par] = dp[par]*10
    time2, MA2, sp2, p2 = OdeReact(anti,tmul,2*tmul,MA1[-1],dp,par)
    
    
    # Third step
    
    dp[par] = dp[par]*2
    time3, MA3, sp3, p3 = OdeReact(anti,2*tmul,3*tmul,MA2[-1],dp,par)
    
    # Fourth step
    
    dp[par] = dp[par]*1.5
    time4, MA4, sp4, p4 = OdeReact(anti,3*tmul,4*tmul,MA3[-1],dp,par)
    
    # Pack all steps into 1 vector / matrix 
    
    MAf = np.row_stack([MA1,MA2,MA3,MA4])
    tf = np.append(np.append(np.append(time1,time2),time3),time4)
    spf = sp1+sp2+sp3+sp4
    pf = p1 + p2 + p3 + p4
    
    #plt.plot(tf,MAf[:,-3:-2]) # plot x3
    #plt.ylabel('Y(t)')
    #plt.xlabel('Time (s)')
    #plt.plot(tf,spf,'--k', linewidth=1)
    #plt.plot(tf,pf,':r', linewidth=2)
    #plt.legend(['x3','Mu',par],loc='best')
    ##plt.title('Near Perfect Adaptation')
    
    plt.plot(tf,MAf[:,-3:]) # plot v and x3
    plt.ylabel('Y(t)')
    plt.xlabel('Time (s)')
    plt.plot(tf,spf,'--k', linewidth=1)
    plt.plot(tf,pf,':r', linewidth=2)
    plt.legend(['x3','Z1','Z2','Mu',par],loc='upper right')
    plt.title('Debugging')
    
    #plt.plot(tf,MAf) # Plot everything
    #plt.ylabel('Y(t)')
    #plt.xlabel('Time (s)')
    #plt.plot(tf,spf,'--k', linewidth=1)
    #plt.legend(['x1','x2','x3','Z1','Z2','Mu'],loc='upper right')
    return 0


def twosteps(contrl):
    # first step
    
    X0 = np.array([0.1, 1, .1, .1, .1])
    time1, MA1, sp1, p1 = OdeReact(contrl,0,tmul,X0,dp,par)
    
    
    # second step
    
    dp[par] = dp[par]* parscale
    time2, MA2, sp2, p2 = OdeReact(contrl,tmul,3*tmul,MA1[-1],dp,par)
    
    # Pack all steps into 1 vector / matrix 
    
    MAf = np.row_stack([MA1,MA2])
    tf = np.append(time1,time2)
    spf = sp1+sp2
    pf = p1 + p2 
    
#    plt.plot(tf,MAf[:,-3:-2]) # plot x3
#    plt.ylabel('Y(t)')
#    plt.xlabel('Time (s)')
    
#    plt.plot(tf,spf,'--k', linewidth=1)
#    plt.plot(tf,pf,':r', linewidth=2)
#    plt.legend(['x3','Mu',par],loc='best')
#    plt.title('System perturbation : Antithetical network')
    #plt.title('Near Perfect Adaptation')
    
    plt.plot(tf,MAf[:,-3:]) # plot v and x3
    plt.ylabel('Y(t)')
    plt.xlabel('Time (s)')
    plt.plot(tf,spf,'--k', linewidth=1)
    plt.plot(tf,pf,':r', linewidth=2)
    plt.legend(['x3','Z1','Z2','Mu',par],loc='upper right')
    plt.title('Debugging')
    
    
    dp[par] = dp[par]/ parscale # return the parameter to original value
    return [tf,spf,pf]

twosteps(anti)
## Hill step
#    
#twosteps(hill)


# Anti step

#twosteps(anti)
#
## intgl step
#
#[tf, spf,pf] = twosteps(intgl)
#
## Pack all steps into 1 vector / matrix 
#
##MAf = np.row_stack([MA1,MA2])
##tf = np.append(time1,time2)
##spf = sp1+sp2
##pf = p1 + p2 
#
#
#plt.plot(tf,spf,'--k', linewidth=1)
#plt.plot(tf,pf,':r', linewidth=2)
#plt.legend(['Antithetical','Integral','Mu',par],loc='best')
#plt.title('Comparision of controller architectures')