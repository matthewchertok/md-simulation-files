# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 12:46:37 2025

@author: zh509
"""

import numpy as np
import matplotlib.pyplot as plt

lam = np.array([0.009219683,0.047941372,0.115048663,0.206341023,0.316084251,0.437383296,0.562616704,0.683915749,0.793658977,0.884951337,0.952058628,0.990780317])
wei = np.array([0.023587668,0.053469663,0.080039164,0.101583713,0.116746268,0.124573523,0.124573523,0.116746268,0.101583713,0.080039164,0.053469663,0.023587668])

dU0 = np.zeros(12) # (U(lambda+d) - U(lambda-d))/2d
dUp = np.zeros(12) # (U(lambda+d) - U(lambda))  /d
dUm = np.zeros(12) # (U(lambda)   - U(lambda-d))/d

for i in range(12):
    output = open(f"500/{i+1}/{i+1}.out")
    lines = output.readlines()
    output.close()
    
    temp_U = [[],[],[]]
    
    iblock = 0
    iline = 0
    while iline < len(lines):
        if lines[iline].startswith('Step'):
            while(1):
                iline = iline + 1
                if lines[iline].startswith('WARNING'):
                    pass
                elif lines[iline].startswith('Loop time of'):
                    break
                else:
                    line = lines[iline].split()
                    temp_U[iblock].append(float(line[2]))
            iblock += 1
        iline += 1
    
    temp_U = np.array(temp_U)
    
    dU0[i] = np.average(temp_U[2]-temp_U[1])/2/0.002
    dUp[i] = np.average(temp_U[2]-temp_U[0])/0.002
    dUm[i] = np.average(temp_U[0]-temp_U[1])/0.002
    
plt.figure(figsize=(4,3))
plt.plot(lam,dU0,'.-',label=r'$\left\langle\frac{U(\lambda+\delta)-U(\lambda-\delta)}{2\delta}\right\rangle$')
plt.plot(lam,dUp,'.-',label=r'$\left\langle\frac{U(\lambda+\delta)-U(\lambda)}{\delta}\right\rangle$')
plt.plot(lam,dUm,'.-',label=r'$\left\langle\frac{U(\lambda)-U(\lambda-\delta)}{\delta}\right\rangle$')
plt.hlines(0,0,1,color='black',linestyles='dashed')
plt.xlabel(r'$\lambda$')
plt.ylabel(r'$\left\langle\frac{dU}{d\lambda}\right\rangle_\lambda$')
plt.xlim(0,1)
plt.legend()
plt.show()

print('Free energy from dU0: ',np.sum(wei*dU0))
print('Free energy from dUp: ',np.sum(wei*dUp))
print('Free energy from dUm: ',np.sum(wei*dUm))
