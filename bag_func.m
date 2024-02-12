function [data_subset,time_subset] = bag_func(data,time,trial_size)
    all_indices = length(data); %number of samples in training data
    num_points = int32( trial_size.* all_indices ); %Number of points in each bag
    subset_indices = sort( randperm(all_indices,num_points) );%random selection of points from 1 to length of data that we then sort from low to high
    data_subset = data([1:size(data,1)],[subset_indices]); %sub-select data points from random indices
    time_subset = time(1,[subset_indices]); %sub-select time points from random indices

end