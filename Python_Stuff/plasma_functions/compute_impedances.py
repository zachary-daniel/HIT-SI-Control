""" This function computes the various impedances of the elements in the HIT-SIU circuits (flux or voltage)"""
from scipy.signal import hilbert, boxcar
import numpy as np
from moving_average import moving_average
import matplotlib.pyplot as plt
def compute_impedances(states, derivs, frequency, window_size = 25, place =None):

    """Parameters:
     states: numpy arrays. This should be organized in the state order of the circuits as follows: i_L1_1, v_cap_1, i_L2_1, etc. 
     derivs: numpy arrays. This isn't really the derivatives, but I'll think of a better name later. These should be ordered as follows: 
     v_L1_1, i_cap_1, v_L2_1
     frequency: double: injector frequency. Will be used to calculate impedance of circuit elements
     window_size: size of the window used for moving average


     returns: 
     L1_arr: array of the impedances of the series inductors in the circuits. We take the capacitor impedance as 96 uF (from the experiment).
     L2_arr:  array of impedances for the flux coils in the flux circuits
     We will include in this imp_L1_1, imp_L2_1 ... . We assume that the mutual coupling between circuits and close neighbors will always be the same, and the weak coupling 
     will always be the same. This list will be of length 8
     M: Double. Strong coupling between circuits
     Mw: Double: strong cpupling between circuits
     R_arr: Since there are no resistive elements in the circuits, the resistive components of the impedance will be used as the resistance 
     in each branch of the circuit.

     """
    if place == None:
        place = int(len(states)/2)

    L1_arr = np.zeros(4) #array of series idnuctor impedances 
    L2_arr = np.zeros(4) #array of flux coil impedances
    R_arr = np.zeros((3,4)) #array of resistors
    count = 0
    #For loop for calculating the individual impedances of all the coils
    for k in range(0,10,3):
        #Need to transpose so the correct dimension is used by hilbert.
        zrat_L1 = hilbert(derivs[:,k].T)/hilbert(states[:,k].T)  #For L1, V_l1/I_L1
        zrat_L2 = hilbert(derivs[:,k+2].T)/hilbert(states[:,k+2].T) #For L2, V_L2,I_L2
        imp_L1 = (moving_average(zrat_L1.imag,window_size) / (2*np.pi*frequency) )[place] # data point is like basically in the middle...
        imp_L2 = (moving_average(zrat_L2.imag,window_size) / (2*np.pi*frequency) )[place]
        R_L1 = moving_average(zrat_L1.real, window_size)[place] #real part of the impedance for L1 branch
        R_L2 = moving_average(zrat_L2.real, window_size)[place] #real part of the impedance for the L2 branch
        #For right now, I'm going to use the L2 resistance as the resistance in the capacitive branch as well
        R_C = R_L2
        R_arr[0,count] = R_L1
        R_arr[1,count] = R_L2
        R_arr[2,count] = R_C
        print(f'imp_L1 = {imp_L1}, imp_L2 = {imp_L2}')
        L1_arr[count] = imp_L1
        L2_arr[count] = imp_L2
        count = count + 1
    
    #Now to calculate the mutual inductances. Coupled equations were computed in Mathematica. As previously stated, we are going to assume that the mutual coupling
    # between all of the circuits is the same, and that the weak coupling is the same as well. we can edit this assumption if need be, but it would be a bit annoying...
    # using the relation that V = L (dI/dt) for an inductor, we can write dI/dT = V_fc/L_fc for each flux coil. 
    # We will compute the mutual coupling, and weak mutual coupling using the first flux circuits
    L21 = L2_arr[0] #Flux coil impedance for first flux circuit
    L22 = L2_arr[1]  #Flux coil impedance for second flux circuit

    #these resistance measurments are sorta hand waved? We can also get these values from the real part of the impedance from the L2 and L1 measurments
     
    R21 = R_arr[1,0] #Ohms
    R22 = R_arr[1,2]  #Ohms. Using the third circuit because the second circuit gives me a negative resistance for some reason...
    R31 = R_arr[2,0] #Ohms
    R32 = R_arr[2,2] #Ohms Using the third circuit because the second circuit gives me a negative resistance for some reason...

    #Need to index all of these at a specific point. For this we choose the sample at place
    x1 = states[place,0]
    x2 = states[place,1]
    x3 = states[place,2]
    x4 = states[place,3]
    x5 = states[place,4]
    x6 = states[place,5]

    # x1 = states[:,0]
    # x2 = states[:,1]
    # x3 = states[:,2]
    # x4 = states[:,3]
    # x5 = states[:,4]
    # x6 = states[:,5]

    # x7 = states[:,6]
    # x8 = states[:,7]
    # x9 = states[:,8]
    # x10 = states[:,9]
    # x11 = states[:,10]
    # x12 = states[:,11]

    x3dot = derivs[:,2]/L2_arr[0] #first flux circuit current derivative
    x6dot = derivs[:,5]/L2_arr[1] #second flux circuit current derivative
    x9dot = derivs[:,8]/L2_arr[2] #third flux circuit current derivative
    x12dot = derivs[:,11]/L2_arr[3] #fourth flux circuit current derivative
    
    # x3dot = derivs[place,2]/L2_arr[0] #first flux circuit current derivative
    # x6dot = derivs[place,5]/L2_arr[1] #second flux circuit current derivative
    # x9dot = derivs[place,8]/L2_arr[2] #third flux circuit current derivative
    # x12dot = derivs[place,11]/L2_arr[3] #fourth flux circuit current derivative


    

    M = -((-R21 * x1 * x12dot - x12dot * x2 + R21 * x12dot * x3 + R31 * x12dot * x3 +  L21 * x12dot * x3dot + R22 * x4 * x9dot + x5 * x9dot - R22 * x6 * x9dot - 
           R32 * x6 * x9dot - L22 * x6dot * x9dot)/(x12dot**2 + x12dot * x6dot - x3dot * x9dot - x9dot**2)) #Strong coupling between nearest circuits
    

    Mw = -((-R21 * x1 * x3dot - x2 * x3dot + R21 * x3 * x3dot + R31 * x3 * x3dot +  L21 * x3dot**2 + R22 * x12dot * x4 + x12dot * x5 - R22 * x12dot * x6 -  R32 * x12dot * x6 
            - L22 * x12dot * x6dot + R22 * x4 * x6dot + x5 * x6dot - R22 * x6 * x6dot - R32 * x6 * x6dot - L22 * x6dot**2 - R21 * x1 * x9dot - x2 * x9dot + R21 * x3 * x9dot + R31 * x3 * x9dot + 
     L21 * x3dot * x9dot)/(-x12dot**2 - x12dot * x6dot + x3dot * x9dot + x9dot**2)) #weak coupling between circuits that are farther apart. 

    plt.figure()
    plt.plot(moving_average(M,window_size))
    plt.figure()
    plt.plot(moving_average(Mw,window_size))
    plt.show

    M_avg = moving_average(M,window_size)
    Mw_avg = moving_average(Mw, window_size)


    

    return L1_arr, L2_arr, M_avg, Mw_avg, R_arr