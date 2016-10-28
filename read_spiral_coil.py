
import numpy as np
import matplotlib.pyplot as plt

file = 'D:\\ykwang\\Data\\SPP simulation\\SPP_list\\Run01\\50R\\kl_1\\99\\999.txt'

data = np.loadtxt(file,dtype='float')
long_chain = data[0:96,0:2]
short_chain = data[96:,0:2]

#PBC_fix
for i in range(1,96):
	for d in range(2):
		diff=long_chain[i,d]-long_chain[i-1,d]
		if abs(diff)>60.:
			if diff>0.:
				long_chain[i,d]-=120.
			else:
				long_chain[i,d]+=120.

plt.plot(long_chain[:,0],long_chain[:,1],'ro')

plt.axis('equal')
plt.show()





