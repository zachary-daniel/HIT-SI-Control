import matplotlib.pyplot as plt
import numpy as np
def make_graph(data,title,legend,impulse,time,log_scale = False): #data should be a list of what one wants plotted,title is a string,
                                  #legend is an array of strings,impulse is a boolean if one is plotting the ringdown of a shot
    '''
    Function for making graphs of circuit data in a particular configuration and order
    data: list of arrays to be plotted on the graph
    title: string to be used as plot title
    legend: array of strings indicating each legend entry. 
    impulse: bool. Changes the x-axis scale if the data being plotted is the tail end of the shot
    time: array. Time vector used when plotting

    '''


    fig, ax = plt.subplots(nrows=3, ncols=4, sharex=True, sharey = 'row', figsize=(12, 10))
    fig.tight_layout()
    
    fig.text(0.5, -0.04, 'Time (ms)', ha='center', fontsize = 20)
    fig.text(-0.04, 0.175, 'Flux Coil Current (Ampere)', va='center', rotation='vertical', fontsize = 20)
    fig.text(-0.04, 0.5, 'Capacitor Voltage (V)', va='center', rotation='vertical', fontsize = 20)
    fig.text(-0.04, 0.825, 'Series Coil Current (Ampere)', va='center', rotation='vertical', fontsize = 20)
    fig.text(.5,1.04,title,ha = 'center',fontsize = 40)
    
    if impulse == False:
        plt.xlim([0,4])


    L1 = [0,3,6,9]
    C = [1,4,7,10]
    L2 = [2,5,8,11]
    
    # colors = ['r','k','--g']
    # alphas = [1,.65]
    for i in range(len(L1)):
        count = 0
        for j in data:
            plt.subplot(3,4,i+1)
            if log_scale == True:
                plt.semilogy(1000*time,np.real(j[:,L1[i]]))
            else:
                plt.plot(1000*time,np.real(j[:,L1[i]]))
            count = count + 1
        plt.grid()


    for i in range(len(C)):
        count = 0
        for j in data:
            plt.subplot(3,4,i+5)
            if log_scale == True:
                plt.semilogy(1000*time,np.real(j[:,C[i]]))
            else:
                plt.plot(1000*time,np.real(j[:,C[i]]))
            count = count + 1
        plt.grid()

    for i in range(len(L2)):
        count = 0
        for j in data:
            plt.subplot(3,4,i+9)
            if log_scale == True:
                plt.semilogy(1000*time,np.real(j[:,L2[i]]))
            else:
                plt.plot(1000*time,np.real(j[:,L1[i]]))

            count = count + 1
        plt.grid()



    fig.legend(legend, fontsize = 15)
    
