function [w_kmeans,e_kmeans,b_avg,atilde,flag] = bop_dmd_func_cluster(num_trials, data,trial_size,rank,imode,train_time,maxiter) %If nothing is passed to trial size, then the default for trial_size = .2

opts = varpro_opts('maxiter', maxiter); %set number of iterations for each run of opt-dmd
data_rows = size(data,1); %number of rows in data. Data should be passed in  as states x samples
w_top_plane = zeros(data_rows*num_trials,ceil(rank/2)); %pre-allocate the space for the eigenvalue  array
if mod(rank,2) == 1 %Determine if there is even or odd rank
    num_points = num_trials*(rank-1)*.5 + num_trials; %number of points in the array
    e_top_plane = zeros(num_points,2); %pre-allocate
else
    e_top_plane = zeros(num_trials*rank*.5,2); %if even. pre-allocate
end

b_total = 0;
batch_size = trial_size; %percentage of data that's used for each run 
shape = size(data);
flag = 'working';
num_points = int16(batch_size.*length(data)); %get the number of points that will be contained in each trial
k = 0;
all_indices = length(data); %length of the data vector
while k < num_trials
    k = k + 1;
    subset_indices = sort( randperm(all_indices,num_points) );%random selection of points from 1 to length of data that we then sort from low to high
    data_subset = zeros([shape(1),num_points]); %pre-allocate an array that's the size of the number of points that will be use in each trial
    time_subset = zeros(size(num_points));
    for v = 1:num_points %for loop to take a random point from our total data matrix and place it into the data_subset
        point = subset_indices(v); %index that's in the random array
        data_subset(:,v) = data(:,point);
        time_subset(v) = train_time(point);
    end
    lbc = [-Inf*ones(rank,1); -Inf*ones(rank,1)]; %Constrain eigenvalues to left half plane
    ubc = [zeros(rank,1); Inf*ones(rank,1)]; 

    copts = varpro_lsqlinopts('lbc',lbc,'ubc',ubc);


    [w,e,b,convergence] = optdmd(data_subset,time_subset,rank,imode,opts,[],[],copts); %call opt-dmd

%Check for convergence
    if convergence == 1 || convergence == 4
        k = k-1; %if convergence fails, redo trial
        

    else
%         w = 0;
%         e = 0;
%         b = 0;


    
        if max(real(e)) > 0 %Check if we get an eigenvalue greater than 0
            flag =  'Not working';
        end
    
    
    
        e_real = real(e); %Break up the real and imaginary parts of eigenvalues
        e_imag = imag(e);
    
        tol = 10; %tolerance for deciding if an eigenvalue is a lone pair
        count = -1; %count for keeping track of where we are in the top plane array
        count_w = 0; %counter for keeping track of row to assign vectors to in w
        start_row_e = ((k-1) * ceil(rank/2)) + 1; %Starting position for the current part of the eigenvalue array 
        start_row_w = (k-1)*data_rows + 1;
        for pos = 1:length(e)
            if e_imag(pos) > 0 || abs(e_imag(pos)) < tol
                count = count + 1;
                count_w = count_w + 1;
                e_top_plane(start_row_e + count,1) = e_real(pos);
                e_top_plane(start_row_e + count,2) = e_imag(pos);
                w_top_plane(start_row_w:start_row_w + data_rows-1,count_w) = w(:,pos); %asign associated eigenvectors to array of eigenvectors in top plane
            end
        end
    end





    
end
%Alright. To preserve that the eigenvalues we get from k-means come in
%complex conjugate pairs, we're only going to do k-means clustering on the
%eigenvalues in the top left half plane (positive imaginary parts). We're
%also going to keep any lone eigenvalues that come out of the opt-dmd
%algorithm when it's run with an odd number of POD modes. We can check this
%as these eigenvalues will have really small imaginary parts, so we just
%threshold them out. The eigenvectors that will be used for averaging can
%be accounted for with the simple fact that the eigenvector for a given
%eigenvalue, is the complex conjugate of the eigenvector for the complex
%conjugate eigenvalue. 

[idx_top_plane,centroids_top_plane] = kmedoids(e_top_plane,ceil(rank/2));
centroids_bottom_plane = zeros(floor(rank/2),2);
w_kmeans = zeros(data_rows,rank);

for k = 1:length(idx_top_plane)
    col = idx_top_plane(k);
    start = floor(k/rank) + 1;
    stop = start + data_rows - 1;
    w_kmeans(:,col) = w_kmeans(:,col) + w_top_plane(start:stop,col);

end

if mod(rank,2) == 1
    [~,I] = min(centroids_top_plane(:,2)); %position of the smallest imaginary part in top plane of centroids. This will be taken to be the lone mode
    centroids_bottom_plane(1:I-1,:) = centroids_top_plane(1:I-1,:);
    centroids_bottom_plane(I:end,:) = centroids_top_plane(I+1:end,:);
    centroids_bottom_plane(:,2) = -centroids_bottom_plane(:,2); %Take conjugate of all of the points in this array

    pos = ceil(rank/2);
    for k = 1:ceil(rank/2)
        if k ~= I
            pos = pos + 1;
            w_kmeans(:,pos) = w_kmeans(:,k);
        end
    end


else
    centroids_bottom_plane = centroids_top_plane;
    centroids_bottom_plane(:,2) = -centroids_bottom_plane(:,2);
    w_kmeans(:,(rank/2) + 1:end) = w_kmeans(:,1:rank/2);
end





centroids = [centroids_top_plane;centroids_bottom_plane];
% w_kmeans = w_kmeans/num_trials;
% start_col = ceil(rank/2)+1;
% 
% end_col = start_col + I-1;







e_kmeans = centroids(:,1) + 1i*centroids(:,2);



b_avg = b_total/num_trials;

atilde = w_kmeans*diag(e_kmeans)*pinv(w_kmeans);
end