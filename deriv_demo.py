
import numpy as np
import matplotlib.pyplot as plt
import math as math 

def smooth_1D(input):
	
	filtered=np.copy(input)
	
	window=np.array([1.,1.,1.])
	window/=sum(window)
	
	for i in range(1,input.size-1):
		filtered[i] = sum(input[i-1:i+2]*window)

	return filtered
	
def label_region_1D(input):
	label = input.copy()
	group = 1 
	for i in range(label.size):
		if label[i] == 1:
			label[i] = group
			if i == (label.size-1):
				break
			if label[i+1] == 0:
				group += 1
	
	return label
	


def flat_detect():
	
	data = np.loadtxt('r_g.txt')
	data_smooth = data.copy()

	continue_threshlod = 25
	smooth_num = 10

	for i in range(smooth_num):
		data_smooth=smooth_1D(data_smooth)
		
	deriv_criteria = np.linspace(0,140,10)

	data_deriv=np.gradient(data_smooth)
	data_flat = np.zeros(data_deriv.size,dtype=int)
	data_flat[np.where(np.abs(data_deriv)<0.03)]=1

	flat_label = label_region_1D(data_flat)
	hist_label = np.histogram(flat_label,bins=np.max(flat_label)+1)[0]
	hist_label = hist_label[1:]

	#from tuple(array) to int, so I need two times of []
	label_max = np.where(hist_label==np.max(hist_label))[0][0]+1
	max_flat_ragion = np.where(flat_label==label_max)[0]
	
	print(max_flat_ragion.size)
	
	if max_flat_ragion.size <= continue_threshlod:
		return 0
	return max_flat_ragion
"""
	plt.subplot(121)
	plt.plot(data,'black')
	plt.plot(max_flat_ragion,data[max_flat_ragion],'r',linewidth=2)

	plt.title('Signal')

	plt.subplot(122)
	plt.plot(data_deriv,'r')
	plt.plot(deriv_criteria,-0.03*np.ones(10))
	plt.title('First Order Deriv')

	plt.show()
"""
	
	

