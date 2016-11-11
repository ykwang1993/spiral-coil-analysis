
import numpy as np
import matplotlib.pyplot as plt
import math as math 
import scipy.signal as signal

def PBC_fix(long_chain):
	for i in range(1,96):
		for d in range(2):
			diff=long_chain[i,d]-long_chain[i-1,d]
			if abs(diff)>60.:
				if diff>0.:
					long_chain[i,d]-=120.
				else:
					long_chain[i,d]+=120.
	return long_chain
	
def persistence_length(long_chain,corr_mode):
	
	vec = long_chain[1:,:]-long_chain[:-1,:]
	vec_length = np.linalg.norm(vec,axis=1)
	vec_length = np.vstack((vec_length,vec_length)).transpose()
	vec_norm = vec/vec_length
	p_l=np.correlate(vec_norm[:,0],vec_norm[:,0],corr_mode) + np.correlate(vec_norm[:,1],vec_norm[:,1],corr_mode)
	p_l=p_l[:int(len(p_l)/2)]
	p_l=p_l[::-1]
	p_l/=np.max(p_l)
	
	return p_l
	
def persistence_length_mod(long_chain):

	vec = long_chain[1:,:]-long_chain[:-1,:]
	vec_length = np.linalg.norm(vec,axis=1)
	vec_length = np.vstack((vec_length,vec_length)).transpose()
	vec_norm = vec/vec_length
	
	vec_0 = np.ones((vec_norm.shape[0],2))*vec_norm[0,:]
	p_l=(vec_0*vec_norm).sum(axis=1)

	p_l/=np.max(p_l)
	
	return p_l
	
def smooth(long_chain):
	
	filtered=np.copy(long_chain)
	
	window=np.array([1.,1.,1.])
	window/=sum(window)
	
	for i in range(1,95):
		filtered[i,0] = sum(long_chain[i-1:i+2,0]*window)
		filtered[i,1] = sum(long_chain[i-1:i+2,1]*window)

	return filtered

	

	
	
def plot_analysis():
	
	outline = 240
	
	plt.subplot(outline+1)
	plt.plot(long_chain[:,0],long_chain[:,1],'ro-')	
	plt.axis('equal')
	plt.title('Zoom In')
	plt.xlabel('x')
	plt.ylabel('y')

	
	plt.subplot(outline+2)
	plt.plot(p_l_f,'b')
	#plt.plot(np.log(p_l_f),'b')
	plt.title('P_l full')
	
	plt.subplot(outline+3)
	plt.plot(p_l_s,'b')
	#plt.semilogy(p_l_s)
	#plt.plot(np.log(p_l_s),'b')
	plt.title('P_l same')
	
	plt.subplot(outline+4)
	plt.plot(p_l_mod,'b')
	plt.title('P_l mod')
	
	plt.subplot(outline+5)
	plt.semilogy(p_l_mod,'b')
	#plt.plot(np.log(p_l_mod),'b')
	plt.title('P_l mod ylog')
	
	
	
	plt.subplot(outline+6)
	
	fit_index = np.array(np.where(p_l_mod>0.)).reshape(-1)
	fit_index = fit_index[1:]-fit_index[:-1] #difference between two index 
	
	print(np.array(np.where(fit_index!=1)[0]).reshape(-1))
	
	#To avoid all the p_l_mod is positive 
	fit_cut_point = np.array(np.where(fit_index!=1)[0]).reshape(-1)
	
	if fit_cut_point.size==0: 
		fit_cut_point=len(p_l_mod)
	else: 
		fit_cut_point=fit_cut_point[0]
	
	#The index means that where the data start to become discontinuous
	fit_data = np.log(p_l_mod[:fit_cut_point+1])
	fit_x = np.arange(len(fit_data))
	A = np.vstack([fit_x, np.ones(len(fit_x))]).T
	fit_result = np.linalg.lstsq(A, fit_data)[0]
	
	plt.plot(fit_data,'b') 
	plt.title('P_l mod ylog fit')
	plt.plot(fit_x,fit_result[0]*fit_x+fit_result[1],'r')
	
	print(fit_result)
	



	plt.show()	




index = '{:03d}'.format(999)
file = 'D:\\ykwang\\Data\\SPP simulation\\SPP_list\\Run01\\15R\\kl_100\\4\\'+index+'.txt'
	
data = np.loadtxt(file,dtype='float')
long_chain = data[0:96,0:2]
short_chain = data[96:,0:2]

long_chain=PBC_fix(long_chain)


for s_i in range(10):
	long_chain=smooth(long_chain)












#persistence length
p_l_f=persistence_length(long_chain,'full')
p_l_s=persistence_length(long_chain,'same')
p_l_mod=persistence_length_mod(long_chain)

plot_analysis()






