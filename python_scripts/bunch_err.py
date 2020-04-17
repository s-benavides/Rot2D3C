from numpy import *
import matplotlib.pyplot as plt

def bunch_err(data,err_ind = -1):
	avg = mean(data)
	err= std(data)/sqrt(len(data))
        avg_list = [avg]
	err_list = [err]
        list = ndarray.tolist(data)[:]
        iter = int(log(len(list))/log(2))
        for w in range(iter):
        	new_list=[]
                while len(list)>1:
                	x=list.pop()
                        y=list.pop()
                        new_list.append((x+y)/2.)
                list=new_list[:]   # the "[:]" is essential list is a copy of new_list.
                avg_list.append(mean(list))
                err_list.append(sqrt((sum((s - avg_list[w+1])**2 for s in list))/((float(len(list)))**2)))
	if err_ind<0:
	        plt.figure()
        	plt.xlabel('# of iterations')
        	plt.ylabel(r"Error")
        	plt.plot(err_list,'.k')
        	plt.show()
		error = input('What is the error? :')
	else:
		error = err_list[err_ind]

	return error
