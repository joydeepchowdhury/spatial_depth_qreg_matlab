This folder contains five MATLAB functions, seven MATLAB programs and four MATLAB datasets. To use them, put all of these in your MATLAB directory, and then do as written in DESCRIPTION OF MATLAB PROGRAMS below.

### The MATLAB functions:
1. crossvalidation.m
2. kernelweights.m
3. spatialquantile.m
4. spatialquantileconfidenceset.m
5. wsdrank.m

### The MATLAB programs calling the above functions:
1. confsetcigar.m
2. confsetsimheter.m
3. confsetsumhes.m
4. qr1.m
5. qr2.m
6. qr3.m
7. spatialdepthdemo.m

### The MATLAB datasets used by the above programs:
1. qr1dataheter.mat
2. qr1datahom.mat
3. Cigar.mat
4. SumHes.mat

## DESCRIPTION OF MATLAB FUNCTIONS:

1. crossvalidation.m

	optimum_h_or_neighborhood_size = crossvalidation(t_vector_X, X_static, t_vector_Y, Y_static, method_for_h, type, Kernel)
	
		Arguments: t_vector_X, X_static, t_vector_Y, Y_static, method_for_h, type, Kernel
			t_vector_X		: 1-by-p vector
			X_static		: n-by-p matrix
			t_vector_Y		: 1-by-q vector
			Y_static		: n-by-q matrix
			method_for_h		: 1 or 2
			type			: string, one of 'pointwise_mean', 'pointwise_median' and 'spatial_median'
			Kernel			: function handle for the underlying kernel
		Output: optimum_h_or_neighborhood_size
			
	**Description:**
	This function calculates the optimum bandwidth or neighborhood size using cross validation. If method_for_h = 1, it returns the optimum bandwidth, else it returns the optimum neighborhood size. X_static is the n-by-p data matrix of n covariate observations, Y_static is the n-by-q data matrix of the corresponding response observations. Three types of estimators can be used for cross validations, viz., the weighted pointwise mean, the weighted pointwise median and the weighted spatial median. To use the weighted pointwise mean, put type = 'pointwise_mean', to use the weighted pointwise median, put type = 'pointwise_median', and to use the weighted spatial median, put type = 'spatial_median'. The numbers n, p and q must be greater than 1. Kernel is a function handle representing the kernel function used to calculate the weights. t_vector_X and t_vector_Y are the respective grids on which the observations in X_static and Y_static are recorded.

2. kernelweights.m

	Weights = kernelweights(x, X_static, t_vector, h, Kernel)
	
		Arguments: x, X_static, t_vector, h, Kernel
			x			: 1-by-p vector
			X_static		: n-by-p matrix
			t_vector		: 1-by-p vector of increasing distinct real numbers
			h			: positive real number
			Kernel			: function handle for the underlying kernel
		Output: Weights
			Weights			: n-by-1 vector
			
	**Description:**
	This function computes the kernel weights for the observation present in the data matrix X_data according to the kernel function Kernel and bandwidth h. X_static is a n-by-p data matrix consisting of n functional observations recorded at p points specified in the vector t_vector. The variable x is a row vector of length p, whose elements are the recorded values of the underlying function at the points in t_vector. The argument named Kernel is a function handle specifying the kernel function, and h is a positive number which is the bandwidth of the kernel function.

3. spatialquantile.m

	Quantile = spatialquantile(Data_original, Weights, u_index, c, t_vector)
	
		Arguments: Data_original, Weights, u_index, c, t_vector
			Data_original		: n-by-p matrix
			Weights			: n-by-1 vector
			u_index			: non-negative integer
			c			: real number in (-1, 1)
			t_vector		: 1-by-p vector of increasing distinct real numbers
		Output: Quantile
			Quantile		: 1-by-p vector
			
	**Description:**
	This function computes the weighted spatial quantile from the data matrix X_data_original with corresponding weights in the variable Weights. The data matrix X_data_original is of dimension n-by-p and contains n functional observations which are recorded on a grid of length p, recorded in the row vector t_vector. The variable Weights is a column vector of length n, whose entries are the weights of the corresponding observations in X_data_original. The variables u_index and c are a positive integer and a real number in (-1, 1) respectively, which together determine the variable u for the spatial quantile computation. For example, if u_index = 1 and c = 0.5, then we compute the weighted spatial quantile corresponding to u equalling to 0.5 times the weighted first principal component of the data in X_static_original.

4. spatialquantileconfidenceset.m

	ConfidenceSet = spatialquantileconfidenceset(Data_original, Weights, u_index, c, t_vector, alpha)

		Arguments: Data_original, Weights, u_index, c, t_vector, alpha
			Data_original		: n-by-p matrix
			Weights			: n-by-1 vector
			u_index			: non-negative integer
			c			: real number in (-1, 1)
			t_vector		: 1-by-p vector of increasing distinct real numbers
			alpha			: (1 - alpha) is the asymptotic coverage probability of the confidence set
		Output: ConfidenceSet
			ConfidenceSet		: 2-by-p matrix of the upper boundary and the lower boundary of the confidence set

	**Description:**
	Computes the asymptotic confidence sets of a spatial quantile.
			
5. wsdrank.m

	rankings = wsdrank(X_to_rank, X_data, X_data_weights, t_vector)
	
		Arguments: X_to_rank, X_data, X_data_weights, t_vector
			X_to_rank		: m-by-p matrix
			X_data			: n-by-p matrix
			X_data_weights	: n-by-1 vector
			t_vector		: 1-by-p vector of increasing distinct real numbers
		Output: rankings
			rankings		: m-by-1 vector
			
	**Description:**
	This function calculates the ranks the observations in X_to_rank according to the decreasing value of the weighted spatial depth with respect to  the data matrix X_data and corresponding weights in the column vector X_data_weights. The dimension of the data matrix X_data is n-by-p, where n is the sample size and p is the dimension of an observation vector. The dimension of the matrix X_to_rank is m-by-p, where m is the number of observations to rank. X_data_weights is a n-by-1 column vector, and t_vector is a row vector recording the grid on which the values of the functional observations underlying X_data and X_to_rank are recorded. The number p must be greater than 1.
			

			

## DESCRIPTION OF MATLAB PROGRAMS:

1. confsetcigar.m

	Produces the plots of the asymptotic confidence sets in the Cigar Data.

2. confsetsimheter.m

	Produces the plots of the asymptotic confidence sets in the heteroscedastic model.

3. confsetsumhes.m

	Produces the plots of the asymptotic confidence sets in the Penn Table Data.

4. qr1.m

	This program generates the plots of the conditional spatial quantiles and the conditional spread measures in the heteroscedastic model considered in the paper.

5. qr2.m

	This program generates the plots of the conditional spatial quantiles and the conditional spread measures in the Cigar Data considered in the paper.

6. qr3.m

	This program generates the plots of the conditional spatial quantiles and the conditional spread measures in the Penn Table considered in the paper.

7. spatialdepthdemo.m

	The program generetes a figure to pictorially demonstrate the concept of the spatial depth.




## DESCRIPTION OF MATLAB DATASETS:

1. qr1dataheter.mat

	This dataset contains a simulated sample from the heteroscedastic model described in the paper.

2. qr1datahom.mat

	This dataset contains a simulated sample from a homoscedastic model.

3. Cigar.mat

	This dataset contains the Cigar Data described in the paper.

4. SumHes.mat

	This dataset contains the Penn Table Data described in the paper.
