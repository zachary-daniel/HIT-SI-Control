from scipy.signal import lsim, dlsim

def compare_operators(operator_list,sys_matrices,input,time,sample_time = None,initial_condition = None):
    """
    operator_list: list of the different matrices that the user would like to compare
    sys_matrices: tuple of B,C,D matrix for each simulation
    input: input signal to be used for each simulation
    time: time vector for the simulation

    returns: trajectory_list: list of each trajectory computed during the simulation
    """

    trajectory_list = []
    residual_list = []
    B = sys_matrices[0]
    C = sys_matrices[1]
    D = sys_matrices[2]
    if sample_time:
        for operator in operator_list:

            sysd = (operator,B,C,D,sample_time)
            trajectory = dlsim(sysd,input.T,time,initial_condition)[-1]
            trajectory_list.append(trajectory)
        return trajectory_list
    else:
        for operator in operator_list:
            sysc = (operator,B,C,D)
            trajectory = lsim(sysc,input.T,time,initial_condition)[-1]
            trajectory_list.append(trajectory)
        return trajectory_list
        
        