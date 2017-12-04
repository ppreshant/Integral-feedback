# -*- coding: utf-8 -*-
"""
Created on Mon Dec 04 02:19:50 2017

@author: user
"""

import matplotlib.pyplot as plt

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

#plt.plot(tf,MAf[:,-3:]) # plot v and x3
#plt.ylabel('Y(t)')
#plt.xlabel('Time (s)')
#plt.plot(tf,spf,'--k', linewidth=1)
#plt.plot(tf,pf,':r', linewidth=2)
#plt.legend(['x3','Z1','Z2','Mu',par],loc='upper right')
#plt.ylim(0,20)
#
#plt.title('Debugging')


plt.plot(tf,MAf) # Plot everything
plt.ylabel('Y(t)')
plt.xlabel('Time (s)')
plt.plot(tf,spf,'--k', linewidth=1)
plt.legend(['x1','x2','x3','Z1','Z2','Mu'],loc='upper right')
plt.ylim(0,20)

plt.title('Debugging')