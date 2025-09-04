function H = generate_dirac_comb(n, m, height)
	% H = generate_dirac_comb(n, m, height)
	% Discretization of -d^2/dx^2 on [0,n] with periodic BCs
	% and a Dirac-comb potential
	%	n: interval is [0, n]
	%	m: discretization points on each interval of unit length
	%	height: height of each spike

	if nargin < 3
		height  = m^2; % spike height
	end    

	N = n*m;       % total points
	h = 1/m;       % stepsize
	
	% second derivative 
	e  = ones(N,1);
	D2 = spdiags([e -2*e e], -1:1, N, N);
	D2(1,end) = 1;
	D2(end,1) = 1;
	L = -(1/h^2) * D2;
	
	% Dirac-comb potential
	V = zeros(N,1);
	for k = 1:n-1
		idx = k*m + 1;
		V(idx) = V(idx) + height;
	end
	
	% hamiltonian
	H = L + spdiags(V,0,N,N);
	end
	