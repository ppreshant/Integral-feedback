# -*- coding: utf-8 -*-
"""
Created on Wed Dec 06 03:00:44 2017

@author: Prashant
"""
import numpy as np
#from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy.linalg as nl

# base parameters
# Rate constants : System parameters
du = 1.0; k2 = 1.0; k3 = 1.0;
g1 = 2.0; g2 = 3.0; g3 = 1.0;

# control parameters
k = .1; n = 1.0; mu = 3.0; theta = 10.0; # general
k4 = 0; k5 = 0; g4 = 0; g5 = 0; # synth and degradation of control species
km = 1; h = 1;  # hill fn
# pack parameters into dictionary
dp = {'k': k,'k2': k2,'k3': k3,'g1': g1,'g2': g2,'g3': g3,'du' : du, 'mu' : mu,'theta':theta,'n' : n, 'km' : km,'h' : h, 'k4' : k4, 'k5' : k5, 'g4' : g4, 'g5' : g5} # dictionary with all parameters


def anti_sp(*param):
    
#    [k3,k2,te,g3,g2,g1,k,du,mu,n] = param # order here is important : print( str(a)).translate(None,"'")
    [mu, g5, g4, g3, g2, g1, h, k, km, n, k3, k2, k5, k4, theta, du] = param;
    z1ss = g1/k * (k3 + g2)/k2 * g3/k3 * mu - du/k ;
    z2ss = theta/n * mu/z1ss ;
    
    Ab = np.array([[-g1,    0,    0,    k,     0],
               [ k2, -(g2 + k3),  0,    0,     0],
               [ 0,     k3,      -g3,   0,     0],
               [ 0,     0,        0,    -n*z2ss,    -n*z1ss],
               [ 0,     0,       theta, -n*z2ss,     n*z1ss]])
    
    
    [eg, egv] = nl.eig(Ab)
    eg = eg[2:]
    stab = int(max(eg) < 0)  
    mcrit = .5
    # Stab will return 1 if all eigen values are negative
                      
    return [stab, mcrit]

def intgl_triv(*param):
    
    [mu, g5, g4, g3, g2, g1, h, k, km, n, k3, k2, k5, k4, theta, du] = param;
    theta = .04
    x3ss = du/g1 * k2/(k3 + g2) * k3/g3;
    
    Ab = np.array([[-g1,    0,        0,    k,  0],  
               [ k2, -(g2 + k3),  0,    0,  0],
               [ 0,     k3,      -g3,   0,  0],
               [ 0,     0,   0,  theta*(mu - x3ss),    0],
               [0,0,0,0,0]])
#    c = np.array([du,    0,        0,    0, 0]) 
    [eg, egv] = nl.eig(Ab)
    eg = eg[:-1] # ignore last eigen value which is always 0
    stab = int(max(eg) < 0) # ignore last el in eg 
    mcrit = int(mu < x3ss)
    # Stab will return 1 if all eigen values are negative
                      
    return [stab, mcrit]

def intgl_sp(*param):
    
    [mu, g5, g4, g3, g2, g1, h, k, km, n, k3, k2, k5, k4, theta, du] = param;
    theta = .04
    vss = g1/k * (k3 + g2)/k2 * g3/k3 * mu - du/k;
    Ab = np.array([[-g1,    0,        0,    k,  0],  
               [ k2, -(g2 + k3),  0,    0,  0],
               [ 0,     k3,      -g3,   0,  0],
               [ 0,     0,  0  ,  -theta*vss,    0],
               [0,0,0,0,0]])
    [eg, egv] = nl.eig(Ab)
    eg = eg[:-1] # ignore last eigen value which is always 0
    stab = int(max(eg) < 0)  # Stab will return 1 if all eigen values are negative
    mcrit = int(mu - vss)
                      
    return [stab, mcrit]

def testeigs(): # plot eigen values rescaling parameters - for debugging
    par = 'mu' # parameter of interest    
    parscale = 1.0 # how much should the parameter be scaled
    dp[par] = dp[par] * parscale
     
    partuple = tuple(dp.values())
#    print ('Trivial',intgl_triv(*partuple))
    print ('Set point',anti_sp(*partuple)[0])
    return 0 
#testeigs()   

npts = 5 # number of values in the logspace

sysparmu = ['du','k2','k3','g1','g2','g3','mu'];
#sysparmu = ['k2']
for par in sysparmu:
    # go over all parameters and plot the stability plot    
        
    logvals = np.logspace(-2,2,npts);
    stability = [.5 for i in logvals]; 
    mcarr = [0 for i in logvals];
    
    for i in range(npts):
        parscale = logvals[i] # how much should the parameter be scaled
        dp[par] = dp[par] * parscale 
        partuple = tuple(dp.values())
        
        
        [stability[i],mcarr[i]] = anti_sp(*partuple);
#        print (stability[i]);
        dp[par] = dp[par] / parscale # unscale parameter to prevent future confusion
    
#    # seeing if mu criterea explains stability
#    plt.plot(stability, mcarr)
#    plt.ylabel('if mu < X3ss')
#    plt.xlabel('Stability')
#    plt.yticks(range(2))
    
    plt.plot(dp[par] * logvals,stability) # plot stability wrt parameter
    plt.ylabel('Stability')
    plt.xlabel('Parameter')
    plt.xscale('log')
#    plt.yticks(range(2))
    
    #plt.legend(['x3','Z1','Z2','Mu',par],loc='upper right')
plt.title('Stability plot of setpoint solution x3 = mu')
plt.legend(sysparmu,loc='best')