function [w_conjugate,e_conjugate] = conjugate_pairs_func(w,e)
    %e is a vector of complex values that have been ordered by real parts
    %from smallest to largest
    %This process is outlined as follows: we will take neigboring
    % eigenvalues (in terms of their real parts), and average the real and
    % imaginary components of these values to get complex conjugate pairs
    %We also need a process to insure that we average the correct
    %eigenvalues
    e_conjugate = zeros(size(e));
    w_conjugate = zeros(size(w));
    diff = 0; %difference between two eigevalues
    tol = .1; %tolerance for deciding whether or not two eigenvalues are pairs
    k = 1;  %counter
    while k < length(e) %cycle through given eigenvalues
        

        e_real = (real(e(k)) + real(e(k+1)) )/2; %break up and average real parts
        diff = abs(real(e(k)) - real(e(k+1)))/ abs(e_real); %difference between two neighboring eigenvalues divided by the averaged real part. 
        %This gives a percentage difference between the two eigenvalues
        %which we will use to determine if two eigenvalues are "close"
                                                       
        if diff > tol %If this diff is larger than the tol we allow, then this eigenvalue is a lone pair and we will continue without 
        %performing any sort of averaging.
            e_conjugate(k) = e(k);
            w_conjugate(:,k) = w(:,k);
            k = k + 1;% Increment by one when the next value is not a complex conjugate
            continue
        end
        e_imag = (abs(imag(e(k))) + abs(imag(e(k+1))) )/2; %break up and average imaginary parts
        e_conjugate(k) = e_real + 1i.*e_imag; % positive conjugate pair goes in first index (a + bi)
        e_conjugate(k+1) = e_real - 1i.*e_imag; % negative conjugate pair goes in second index (a - bi)
        
        %Now we repeat the same process for the eigenvectors. Two complex
        %conjugate eigenvalues have eigenvectors that are complex
        %conjugates
        w_real = (real(w(:,k)) + real(w(:,k+1)))/2; %break up real and imaginary parts
        w_imag = (abs(imag(w(:,k))) + abs(imag(w(:,k+1))) )/2;
        w_conjugate(:,k) = w_real + 1i.*w_imag; %first position is the positive conjugate pair
        w_conjugate(:,k+1) = w_real - 1i.*w_imag; %second position is the negative conjugate pair

        k = k + 2; %increment by 2 when the next value is a complex conjugate
    end
    %When we exit the while loop, if the last eigenvalue is a lone pair it
    %will not be assigned. We check below if this is the case, and assign a
    %value if necessary
    if k == length(e) && e_conjugate(k) == 0 && e(k) ~= 0
        e_conjugate(k) = e(k);
        w_conjugate(:,k) = w(:,k);
    end
end