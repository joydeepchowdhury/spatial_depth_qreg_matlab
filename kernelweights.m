function Weights = kernelweights(x, X_static, t_vector, h, Kernel)

% This function computes the kernel weights for the observation present in
% the data matrix X_data according to the kernel function Kernel and
% bandwidth h. X_static is a n-by-p data matrix consisting of n functional
% observations recorded at p points specified in the vector t_vector. The
% variable x is a row vector of length p, whose elements are the recorded
% values of the underlying function at the points in t_vector. The argument
% named Kernel is a function handle specifying the kernel function, and h
% is a positive number which is the bandwidth of the kernel function.

Distance_vector = trapz(t_vector, (ones(size(X_static,1),1) * x - X_static).^2, 2);

Weights = Kernel(Distance_vector / h);
Weights = Weights / sum(Weights);

end