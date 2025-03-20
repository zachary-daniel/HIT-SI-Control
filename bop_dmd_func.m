function [w_avg,e_avg,b_avg,e_var,atilde, atilde_cell_arr] = bop_dmd_func(training_data,time,r,trial_size,num_trials,varargin) %If nothing is passed to trial size, then the default for trial_size = .2
%traning_data: given as states X samples. 
%r: rank of the opt-dmd model
%time: time vector for the training data
%trial_size: percent of training data included in each bag, limited to a
%decimal between 0 and 1.
%num_trials: number of runs of opt-dmd 
%varargin: {max_iter, tol, e_init,complex_conjugate}: 
% max_iter: maximum number of iterations of the
%opt_dmd code before stopping. tol: error tol for convergence of opt_dmd,
%e_init: guess for the initial eigenvalues of opt-dmd
%complex_conjugate: bool. If true, will execute conjugate_pairs_func and
%insure that eigenvalues and eigenvectors come in complex conjugate pairs
imode = 2; %Fit data to POD modes
tol = 1e-10;
maxiter = 30;
complex_conjugate = false;

e_init = [];
if nargin > 5 && ~isempty(varargin{1})
    maxiter = varargin{1};
end
if nargin > 5 && ~isempty(varargin{2})
    tol = varargin{2};
end
if nargin > 5 && ~isempty(varargin{3})
    e_init = varargin{3};
end
if nargin > 5 && ~isempty(varargin{4})
    complex_conjugate = varargin{4};
end

opts = varpro_opts('maxiter', maxiter,'tol',tol);
num_states = size(training_data,1); %number of states in training data
%Following 3 variables will all be of complex type.
w_sum = zeros(num_states,r); %sum of the eigenvectors from each opt-dmd run
e_sum = zeros(r,1); %sum of eigenvalues from each opt-dmd run
b_sum = zeros(r,1); %sum of weights from each opt-dmd run
e_arr_real = zeros(r,num_trials);
e_arr_imag = zeros(r,num_trials);
e_var = zeros(r,2);
successful_trials = 0; %number of successful opt-dmd runs.
atilde_cell_arr = cell(1,num_trials);
while successful_trials < num_trials
    [data_subset,time_subset] = bag_func(training_data,time,trial_size); %get random data and time points to be used in each run of bop-dmd
    
    lbc = [-Inf*ones(r,1); -Inf*ones(r,1)]; %Constrain eigenvalues to left half plane
    ubc = [zeros(r,1); Inf*ones(r,1)]; 
    copts = varpro_lsqlinopts('lbc',lbc,'ubc',ubc);
    [w,e,b,convergence,a_opt] = optdmd(data_subset,time_subset,r,imode,opts,e_init,[],copts);
    e_arr_real(:,successful_trials+1) = real(e);
    e_arr_imag(:,successful_trials+1) = imag(e);
    
    %Initially, I evaluated if the opt-dmd algorithm had converged, but
    %honestly, it never converges for this problem anyways. Removing this
    %check has only improved results. Key is to do an initial run of
    %opt-dmd, and then let the bagging handle just lightly shifting the
    %eigenvalues and vectors. This improves the overall error
    if  convergence ~= 4 
        [e_sort,sorted_indices] = sort(e,'ComparisonMethod','real'); %Sort eigenvalues by increasing real part.
        %add to the running total of all e,b,and w matrices
        e_sum = e_sum + e_sort;
        w_sum = w_sum + w(:,sorted_indices);
        b_sum = b_sum + b(sorted_indices);
        successful_trials = successful_trials + 1;
        atilde_cell_arr{successful_trials} = w*diag(e)*pinv(w);
    else 
       continue
    end
end
e_var(:,1) = var(e_arr_real,0,2); %Take variance of real part of eigs
e_var(:,2) = var(e_arr_imag,0,2); %Variance of imag part of eigs
%below function will ensure that the returned eigenvalues come in complex
%conjugate pairs. 
if complex_conjugate == true
    [w_sum_conjugate,e_sum_conjugate] = conjugate_pairs_func(w_sum,e_sum);
    w_avg = w_sum_conjugate/num_trials;
    e_avg = e_sum_conjugate/num_trials;
    b_avg = b_sum/num_trials;
else
    w_avg = w_sum/num_trials;
    e_avg = e_sum/num_trials;
    b_avg = b_sum/num_trials;
end
atilde = w_avg*diag(e_avg)*pinv(w_avg);

end