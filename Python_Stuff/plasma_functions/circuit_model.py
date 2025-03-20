''' This function will be used to load vacuum models of the HIT-SIU flux circuits with different circuit parameters. 
This will return a single dynamics matrix, A '''
import numpy as np

def circuit_model(L1_arr, L2_arr, Cap_arr, R_arr, M_arr, Vacuum = True):
    ''' 
    L1_arr: array of inductances corresponding to each series inductor from 1-4
    L2_arr: array of inductances corresponding to each flux coil from 1-4. Pass in a fifth entry for the plasma inductance if desired
    Cap_arr: array of capacitances corresponding to each capacitor from flux circuits 1-4
    R_arr: 2d array of resitances for each circuit. will be of dim 3 x 4. 3 resistances for each circuit, R1, R2, R3, and 4 total circuits.
    If a plasma is present, add another row to this array so it is 4x4, and add in the plasma resistance as the last column
    M_arr: array of mutual inductances. Pass in a thrid value for plassma mutual if desired
    '''
    # Will use the relation for a linear system A*xdot = B* x => xdot = A^-1 * B x to get our vacuum model



    # coeff in front of the states x in vacuum
    state_coeff = np.array([
        [-(R_arr[0,0]+R_arr[1,0])/L1_arr[0], -1/L1_arr[0], R_arr[1,0]/L1_arr[0], 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1/Cap_arr[0], 0, -1/Cap_arr[0], 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-R_arr[1,0], -1, R_arr[2,0], 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, (-1/L1_arr[1])*(R_arr[0,1]+R_arr[1,1]), -1/L1_arr[1], R_arr[1,1]/L1_arr[1], 0, 0, 0, 0, 0,0],
        [0, 0, 0, 1/Cap_arr[1], 0, -1/Cap_arr[1], 0, 0, 0, 0, 0, 0],
        [0, 0, 0, -R_arr[1,1], -1, R_arr[2,1], 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, (-1/L1_arr[2])*(R_arr[0,2]+R_arr[1,2]), -1/L1_arr[2], R_arr[1,2]/L1_arr[2], 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1/Cap_arr[2], 0, -1/Cap_arr[2], 0, 0, 0],
        [0, 0, 0, 0, 0, 0, -R_arr[1,2], -1, R_arr[2,2], 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, (-1/L1_arr[3])*(R_arr[0,3]+R_arr[1,3]), -1/L1_arr[3], R_arr[1,3]/L1_arr[3]],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1/Cap_arr[3], 0, -1/Cap_arr[3]],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, -R_arr[1,3], -1, R_arr[2,3]]
    ])

    # coeff in front of state derivatives x_dot in vacuum
    state_derivative_coeff = np.array([[1,0,0,0,0,0,0,0,0,0,0,0],
                            [0,1,0,0,0,0,0,0,0,0,0,0],
                            [0,0,-L2_arr[0],0,0,-M_arr[0],0,0,-M_arr[1],0,0,-M_arr[0]],
                            [0,0,0,1,0,0,0,0,0,0,0,0],
                            [0,0,0,0,1,0,0,0,0,0,0,0],
                            [0,0,-M_arr[0],0,0,-L2_arr[1],0,0,-M_arr[0],0,0,-M_arr[1]],
                            [0,0,0,0,0,0,1,0,0,0,0,0],
                            [0,0,0,0,0,0,0,1,0,0,0,0],
                            [0,0,-M_arr[1],0,0,-M_arr[0],0,0,-L2_arr[2],0,0,-M_arr[0]],
                            [0,0,0,0,0,0,0,0,0,1,0,0],
                            [0,0,0,0,0,0,0,0,0,0,1,0],
                            [0,0,-M_arr[0],0,0,-M_arr[1],0,0,-M_arr[0],0,0,-L2_arr[3]]],
                            )
    #B in vacuum case 
    B_analytic = np.array([[1/L1_arr[0], 0, 0, 0],
                    [0,0,0,0],
                    [0,0,0,0],
                    [0,1/L1_arr[1],0,0],
                    [0,0,0,0],
                    [0,0,0,0],
                    [0,0,1/L1_arr[2],0],
                    [0,0,0,0],
                    [0,0,0,0],
                    [0,0,0,1/L1_arr[3]],
                    [0,0,0,0],
                    [0,0,0,0]]
                    )
    #C in vacuum case
    C_analytic = np.array(
    [[0,0,1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]]
    )
    if Vacuum == False:

        # coeff in front of the states x plasma case
        state_coeff = np.array([
            [-(R_arr[0,0]+R_arr[1,0])/L1_arr[0], -1/L1_arr[0], R_arr[1,0]/L1_arr[0], 0, 0, 0, 0, 0, 0, 0, 0, 0,0],
            [1/Cap_arr[0], 0, -1/Cap_arr[0], 0, 0, 0, 0, 0, 0, 0, 0, 0,0],
            [-R_arr[1,0], -1, R_arr[2,0], 0, 0, 0, 0, 0, 0, 0, 0, 0,0],
            [0, 0, 0, (-1/L1_arr[1])*(R_arr[0,1]+R_arr[1,1]), -1/L1_arr[1], R_arr[1,1]/L1_arr[1], 0, 0, 0, 0, 0,0,0],
            [0, 0, 0, 1/Cap_arr[1], 0, -1/Cap_arr[1], 0, 0, 0, 0, 0, 0,0],
            [0, 0, 0, -R_arr[1,1], -1, R_arr[2,1], 0, 0, 0, 0, 0, 0,0],
            [0, 0, 0, 0, 0, 0, (-1/L1_arr[2])*(R_arr[0,2]+R_arr[1,2]), -1/L1_arr[2], R_arr[1,2]/L1_arr[2], 0, 0, 0,0],
            [0, 0, 0, 0, 0, 0, 1/Cap_arr[2], 0, -1/Cap_arr[2], 0, 0, 0,0],
            [0, 0, 0, 0, 0, 0, -R_arr[1,2], -1, R_arr[2,2], 0, 0, 0,0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, (-1/L1_arr[3])*(R_arr[0,3]+R_arr[1,3]), -1/L1_arr[3], R_arr[1,3]/L1_arr[3],0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 1/Cap_arr[3], 0, -1/Cap_arr[3],0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, -R_arr[1,3], -1, R_arr[2,3],0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, R_arr[0,4]]
        ])

        # coeff in front of state derivatives x_dot plasma case
        state_derivative_coeff = np.array([
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, -L2_arr[0], 0, 0, -M_arr[0], 0, 0, -M_arr[1], 0, 0, -M_arr[0], -M_arr[2]],
            [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, -M_arr[0], 0, 0, -L2_arr[1], 0, 0, -M_arr[0], 0, 0, -M_arr[1], -M_arr[2]],
            [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
            [0, 0, -M_arr[1], 0, 0, -M_arr[0], 0, 0, -L2_arr[2], 0, 0, -M_arr[0], -M_arr[2]],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, -M_arr[0], 0, 0, -M_arr[1], 0, 0, -M_arr[0], 0, 0, -L2_arr[3], -M_arr[2]],
            [0, 0, -M_arr[2], 0, 0, -M_arr[2], 0, 0, -M_arr[2], 0, 0, -M_arr[2], -L2_arr[4]]
        ])


        #plasma case
        B_analytic = np.array([[1/L1_arr[0], 0, 0, 0],
                        [0,0,0,0],
                        [0,0,0,0],
                        [0,1/L1_arr[1],0,0],
                        [0,0,0,0],
                        [0,0,0,0],
                        [0,0,1/L1_arr[2],0],
                        [0,0,0,0],
                        [0,0,0,0],
                        [0,0,0,1/L1_arr[3]],
                        [0,0,0,0],
                        [0,0,0,0],
                        [0,0,0,0]]
                        )
        #plasma case
        C_analytic = np.array(
        [[0,0,1, 0, 0, 0, 0, 0, 0, 0, 0, 0,0],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,0],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,0]]
        )

    D_analytic = np.zeros((C_analytic.shape[0], C_analytic.shape[0]))


    A_analytic = np.matmul(np.linalg.inv(state_derivative_coeff), state_coeff)
    return A_analytic, B_analytic, C_analytic, D_analytic