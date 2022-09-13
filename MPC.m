function [nextInput, new_switch_position, results] = MPC(currentStates, referencePoint, sys_d, switch_position, dt, horizon, current_t, Amplitude, J, Runtime)
    if switch_position == 1             %Check which siwtch states are available
        permissible_states = [1,0];
    elseif switch_position == -1
        permissible_states = [-1,0];
    else 
        permissible_states = [-1,0,1];
    end
    if  (current_t+horizon) > Runtime
        t_window = (current_t):dt:Runtime;
    else
        t_window = current_t:dt:(current_t+horizon-dt); %Determine time window
    end
    y = zeros(length(t_window(1,:)), 12);
    for i = 1:length(permissible_states)
        input = ones(size(t_window,2),1) * (permissible_states(i)*(Amplitude));
        divider = length(input);
        input(1:round(divider/8)-1) = 0;
        input(divider-round(divider/8)+2:divider) = 0;
        y(:, ((4*i)-3) :4*i) = lsim(sys_d, [input input input input], t_window, currentStates);
        output = y(end, (i.^2)+1);
        cost(i) = J(referencePoint,output, switch_position, permissible_states(i));
    end
    [M,I] = min(cost);
    results = y(:, ((4*I)-3): 4*I);
    [nextInput] = ones(size(t_window,2),1) * permissible_states(I)*Amplitude; 
    divider = length(nextInput);
     nextInput(1:round(divider/8)-1) = 0;
     nextInput(divider-round(divider/8)+1:divider) = 0;
    new_switch_position = permissible_states(I);
end
