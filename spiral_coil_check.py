
import numpy as np
import matplotlib.pyplot as plt
import math as math 
import deriv

def radius_gyration(long_chain):
	r_g= math.sqrt(((long_chain[:,0]-long_chain[:,0].mean())**2 + (long_chain[:,1]-long_chain[:,1].mean())**2).mean())
	return r_g

def cumulative_angle(long_chain):
	
	vec = long_chain[1:,:]-long_chain[:-1,:]

	a = vec[:-1,:]
	b = vec[1:,:]

	angle = np.arccos((a[:,0]*b[:,0]+a[:,1]*b[:,1])/(np.sqrt(a[:,0]**2+a[:,1]**2)*np.sqrt(b[:,0]**2+b[:,1]**2)))*180/math.pi
    
	#determein the +,- of angle
	cross = a[:,0]*b[:,1]-a[:,1]*b[:,0]
    #If (cross GT 0.0) Then angle[j+1] = -angle[j+1]
	angle[np.where(cross>0)] *= -1 
	
	return angle

	
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
	plt.plot(long_chain[:,0],long_chain[:,1],'ro-')	
	plt.axis('equal')
	plt.title('Zoom In')
	plt.xlabel('x')
	plt.ylabel('y')
	
	plt.subplot(outline+3)
	plt.plot(time_index,L_ee,'b')
	plt.title('End to End Distance')
	plt.xlabel('t')
	plt.ylabel('Length')
	
	plt.subplot(outline+4)
	plt.plot(vec_angle_time,vec_angle_cumulative,'b')
	plt.title('End to End Cumulative Angle')
	plt.xlabel('t')
	plt.ylabel('Degree')
	
	plt.subplot(outline+5)
	plt.plot(time_index,cumu_vec_t,'b')
	plt.title('Cumulative Angle')
	plt.xlabel('t')
	plt.ylabel('Degree')
	
	plt.subplot(outline+6)
	plt.plot(cumu_vec,'b')
	plt.title('Cumulative Angle of last frame')
	plt.xlabel('# of beads')
	plt.ylabel('Degree')

	plt.subplot(outline+7)
	plt.plot(time_index,r_g,'b')	
	plt.title('Radius of Gyration')
	plt.xlabel('t')
	plt.ylabel('Length')

	flat = deriv.flat_detect()
	plt.plot((flat+1)*8-1,r_g[flat],'r',linewidth=2)
		

		

	
	
	"""
	plt.subplot(outline+7)
	#plt.plot(p_l_f,'b')
	plt.plot(np.log(p_l_f),'b')
	plt.title('P_l full')
	
	plt.subplot(outline+8)
	#plt.plot(p_l_s,'b')
	#plt.semilogy(p_l_s)
	plt.plot(np.log(p_l_s),'b')
	plt.title('P_l same')
	"""


	plt.show()	

	

r_g = np.empty([125],dtype='float')
L_ee = np.empty([125],dtype='float')
vec_ee = np.empty([125,2],dtype='float')
cumu_vec_t = np.empty([125],dtype='float')

time_index=np.arange(125)
time_index = (time_index+1)*8-1



for t in range(125):

	index = '{:03d}'.format(time_index[t])
	file = 'D:\\ykwang\\Data\\SPP simulation\\SPP_list\\Run01\\15R\\kl_20\\2\\'+index+'.txt'
	
	
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
	#radius of gyration
	r_g[t]=radius_gyration(long_chain)

	
	#end to end distance
	L_ee[t] = math.sqrt((long_chain[0,0]-long_chain[95,0])**2+(long_chain[0,1]-long_chain[95,1])**2)
	
	
	
	#store the end to end vector
	vec_ee[t,:] = [long_chain[0,0]-long_chain[95,0],long_chain[0,1]-long_chain[95,1]]
	
	#cumulative angle
	cumu_vec=cumulative_angle(long_chain)
	cumu_vec=np.cumsum(cumu_vec)
	cumu_vec_t[t] = cumu_vec[-1]
	
	
	
	print(t)

	

#cumulative angle change in one frame	
cumu_vec=cumulative_angle(long_chain)
cumu_vec=np.cumsum(cumu_vec)
	
#end to end cumulative angle 
a_end = vec_ee[:-2,:]
b_end = vec_ee[1:-1,:]
vec_angle = np.arccos((a_end[:,0]*b_end[:,0]+a_end[:,1]*b_end[:,1])/(np.sqrt(a_end[:,0]**2+a_end[:,1]**2)*np.sqrt(b_end[:,0]**2+b_end[:,1]**2)))*180/math.pi
vec_angle = np.abs(vec_angle)
vec_angle_cumulative = np.cumsum(vec_angle)
vec_angle_cumulative = np.hstack((0,vec_angle_cumulative))
vec_angle_time = np.linspace(0,999,vec_angle_cumulative.shape[0])

#persistence length
p_l_f=persistence_length(long_chain,'full')
p_l_s=persistence_length(long_chain,'same')


plot_analysis()






