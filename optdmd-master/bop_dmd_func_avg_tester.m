%Fuck this shit
function [w_avg,e_avg,b_avg,atilde,flag] = bop_dmd_func_avg(num_trials, data,trial_size,rank,imode,train_time,maxiter,tol,med,e_init) 
opts = varpro_opts('maxiter', maxiter,'tol',tol); %set number of iterations for each run of opt-dmd
data_rows = size(data,1); %number of rows in data. Data should be passed in as states x samples
w_total = zeros(data_rows,rank); %pre-allocate the space for the eigenvector  array
e_med_real = zeros(num_trials,rank); %pre-allocate space for real part of eigenvalue array
e_med_imag = zeros(num_trials,rank);%pre-allocate space for imag part of eigenvalue array

b_total = 0; %weights for opt-dmd
batch_size = trial_size; %percentage of data that's used for each run 
shape = size(data);
flag = 'working'; %flag used to indicate if an unstable eigenvalue gets through
num_points = int16(batch_size.*length(data)); %get the number of points that will be contained in each trial
k = 0; %index used to indicate what trial the bagging loop is on
all_indices = length(data); %length of the data vector
%Begin loop
while k < num_trials
    k = k + 1; 
    subset_indices = sort( randperm(all_indices,num_points) );%random selection of points from 1 to length of data that we then sort from low to high
    data_subset = zeros([shape(1),num_points]); %pre-allocate an array that's the size of the number of points that will be use in each trial
    time_subset = zeros(size(num_points)); %corresponding time points for randomly selected data points
    for v = 1:num_points %for loop to take a random point from our total data matrix and place it into the data_subset
        point = subset_indices(v); %index that's in the random array
        data_subset(:,v) = data(:,point);
        time_subset(v) = train_time(point);
    end
    lbc = [-Inf*ones(rank,1); -Inf*ones(rank,1)]; %Constrain eigenvalues to left half plane
    ubc = [zeros(rank,1); Inf*ones(rank,1)]; 

    copts = varpro_lsqlinopts('lbc',lbc,'ubc',ubc);


    [w,e,b,convergence] = optdmd(data_subset,time_subset,rank,imode,opts,[e_init],[],copts); %call opt-dmd

%Check for convergence
    if convergence == 1 || convergence == 4
        k = k-1; %if convergence fails, redo trial
        

    else %Successful convergence proceeds


    
        if max(real(e)) > 0 %Check if we get an eigenvalue greater than 0
            flag =  'Not working';
        end
    


        [sort_e,I] = sort(real(e));
        w = w(:,I);
        e = e(I);
        e_med_real(k,end) = real(e(end));%Assign lone eigenpair to end of these arrays to preserve these modes
        e_med_imag(k,end) = imag(e(end));
        w_total(:,end) = w(:,end); 
        b_total = b_total + b;

        [~,I] = min(abs(imag(e))); %Find the mode with the slowest oscillation to designate as the lone eigenvalue for future steps
        if mod(rank,2) == 1 && I ~= length(e) %if this eigenvalue is not the last eigenvlaue, go through the steps to reassign it
            temp = e(I);
            e(I) = e(end);
            e(end) = temp;
            temp_w = w(:,I);
            w(:,I) = w(end);
            w(:,end) = temp_w;


        end

        e_real = real(e); %Break up the real and imaginary parts of eigenvalues
        e_imag = imag(e);
        for pos = 1:rank%1:2:(floor(length(e)/2)*2)-1 %cycle through till second to last array element
%             e_real_mean = (e_real(pos) + e_real(pos + 1))/2; %mean the real parts of the eigevnalues at neighboring positions
%             e_imag_mean = (abs(e_imag(pos)) + abs(e_imag(pos+1)))/2; %mean the abs of the imag parts of eigevnalues at neighboring positions
%             e_med_real(k,pos) = e_real_mean; %assign to the median array
%             e_med_real(k,pos+1) = e_real_mean;
%             e_med_imag(k,pos) = e_imag_mean;
%             e_med_imag(k,pos+1) = -e_imag_mean; %assign negative imaginary component to next spot in median array
%             w_real_mean = (real(w(:,pos)) + real(w(:,pos+1)))/2;
%             w_imag_mean = (abs(imag(w(:,pos))) + abs(imag(w(:,pos+1))))/2; %repeat for eigenvectors
%             w_total(:,pos) = w_total(:,pos) + w_real_mean +1i*w_imag_mean; %add eigenvectors to the total
%             w_total(:,pos+1) = w_total(:,pos+1) + w_real_mean -1i*w_imag_mean; %add conjugate eigenvector to the total at the next position

            if abs(e_imag(pos)) < 100
                e_imag(pos) = 0;
                w(:,pos) = real(w(:,pos));
            end
            e_med_real(k,pos) = e_real(pos);
            e_med_imag(k,pos) = e_imag(pos);
            w_total(:,pos) = w_total(:,pos) + w(:,pos);
        end

        
    end


end


w_avg = w_total/(num_trials); %average  eigenvectors


e_avg = (sum(e_med_real,1)/num_trials) + 1i*(sum(e_med_imag,1)/num_trials); 
e_med = median(e_med_real,1) + 1i*(median(e_med_imag,1)); %add up real and imag parts
b_avg = b_total/num_trials;
if med == false
    atilde = (w_avg)*diag(e_avg)*(pinv(w_avg));
else
    atilde = (w_avg)*diag(e_med)*(pinv(w_avg));
end
end
