import numpy as np

from matplotlib import pyplot as plt

from BDTSM import *

Qx=[1,2.8,5.2,9,7.2]
Qy=[2.8,5.2,9,4.8,1.2]

T1x=[1,2,3.4,6.6,8.6]
T1y=[0.8,3,7.2,7.4,1.4]

T2x=[1.4,3,10,5.4,8.6]
T2y=[2.8,7.0,9.4,4.4,2.2]
Q=np.hstack((np.array([Qx]).T,np.array([Qy]).T))
T1=np.hstack((np.array([T1x]).T,np.array([T1y]).T))
T2=np.hstack((np.array([T2x]).T,np.array([T2y]).T))
print('Q:'+str(round(SIM(Q,Q,3),4)))
print('T1:'+str(round(SIM(Q,T1,3),4)))
print('T2:'+str(round(SIM(Q,T2,3),4)))
plt.plot(Qx,Qy,'o-',color='red',label='Q')
plt.plot(T1x,T1y,'o-',color='green',label='T1')
plt.plot(T2x,T2y,'o-',color='blue',label='T2')
plt.legend()
plt.show()