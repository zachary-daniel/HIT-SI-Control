function [nextInput, new_switch_position, results] = MPC(currentStates, referencePoint, sys_d, switch_position, dt, horizon, current_t, Amplitude, J)
    if switch_position == 1             %Check which siwtch states are available
        permissible_states = [1,0];
    elseif switch_position == -1
        permissible_states = [-1,0];
    else 
        permissible_states = [-1,0,1];
    end
    t_window = current_t:dt:current_t+horizon; %Determine time window
    for i = 1:length(permissible_states)
        input = ones(size(t_window,2)) * (permissible_states(i)*(Amplitude));
        input = diag(input);
        y = lsim(sys_d, [input input input input], t_window, currentStates); 
        cost(i) = J(referencePoint,y(end), switch_position, permissible_states(i));
        if cost(i) == min(cost)
            results = y;
        end
    end
    [M,I] = min(cost);
    [nextInput] = ones(size(t_window,2)) * permissible_states(I)*Amplitude; 
    nextInput = diag(nextInput);
    new_switch_position = permissible_states(I);
end
