import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)
cmap="RdBu_r"





Fore_time_perm = np.loadtxt('Forward_flow_v_time.dat')
Back_time_perm = np.loadtxt('Back_flow_v_time.dat')

#### make time series #####

Fore_time_series = [i for i in range(len(Fore_time_perm))]
Back_time_series = [i for i in range(len(Back_time_perm))]

plt.plot(Fore_time_perm/float(100),Fore_time_series,label='Forwards')
plt.plot(Back_time_perm/float(100),Back_time_series,label='Backwards')
plt.legend()
plt.savefig('Permeation_forward_back.png',bbox_inches='tight')
plt.show()

Time_perm_combined = np.sort(np.concatenate((Fore_time_perm,Back_time_perm,[0])))

print(Time_perm_combined)
Combined_time_series = [i for i in range(len(Time_perm_combined))]
plt.plot(Time_perm_combined/float(100),Combined_time_series)
plt.xticks([0,250,500,750,1000])
plt.savefig('Permeation_combined.png',bbox_inches='tight')
plt.show()
