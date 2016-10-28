
import numpy as np
import matplotlib.pyplot as plt
import math as math 


	
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
	
	
	
	
def plot_analysis():
	
	outline = 240
	
	plt.subplot(outline+1)
	plt.plot(short_chain[:,0],short_chain[:,1],'o',markersize=2.5,color='red', markeredgewidth=0.0)
	plt.plot(long_chain[:,0],long_chain[:,1],'o',markersize=2.5,color='black', markeredgewidth=0.0)	
	plt.title('Last Frame')
	plt.xlabel('x')
	plt.ylabel('y')
	plt.xlim(0,120)
	plt.ylim(0,120)
	
	plt.subplot(outline+2)
	plt.plot(long_chain[:,0],long_chain[:,1],'ro')	
	plt.axis('equal')
	plt.title('Zoom In')
	plt.xlabel('x')
	plt.ylabel('y')
	
	

	plt.subplot(outline+3)
	plt.plot(time_index,r_g,'b')	
	plt.title('Radius of Gyration')
	plt.xlabel('t')
	plt.ylabel('Length')

	plt.subplot(outline+4)
	plt.plot(time_index,L_ee,'b')
	plt.title('End to End Distance')
	plt.xlabel('t')
	plt.ylabel('Length')	

	plt.subplot(outline+5)
	plt.plot(vec_angle_time,vec_angle_cumulative,'b')
	plt.title('End to End Cumulative Angle')
	plt.xlabel('t')
	plt.ylabel('Degree')	

	plt.subplot(outline+6)
	plt.plot(cumu_vec,'b')
	plt.title('Cumulative Angle of last frame')
	plt.xlabel('# of beads')
	plt.ylabel('Degree')
	
	plt.subplot(outline+7)
	plt.plot(p_l_f,'b')
	#plt.plot(np.log(p_l_f),'b')
	plt.title('P_l full')
	
	plt.subplot(outline+8)
	plt.plot(p_l_s,'b')
	#plt.semilogy(p_l_s)
	#plt.plot(np.log(p_l_s),'b')
	plt.title('P_l same')



	plt.show()	



for t in range(125):

	index = '{:03d}'.format(t)
	file = 'D:\\ykwang\\Data\\SPP simulation\\SPP_list\\Run01\\00R\\kl_100\\2\\'+index+'.txt'
	
	
	data = np.loadtxt(file,dtype='float')
	long_chain = data[0:96,0:2]
	short_chain = data[96:,0:2]

	#PBC fix
	for i in range(1,96):
		for d in range(2):
			diff=long_chain[i,d]-long_chain[i-1,d]
			if abs(diff)>60.:
				if diff>0.:
					long_chain[i,d]-=120.
				else:
					long_chain[i,d]+=120.
	
	print(t)


#persistence length
p_l_f=persistence_length(long_chain,'full')
p_l_s=persistence_length(long_chain,'same')


plot_analysis()






