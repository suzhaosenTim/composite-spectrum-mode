function polariz = asymptot_high_contrast(n,rad, epsilon1, delta)

% function of computing polarizability under horizontal incident field
% for a pair of disks. (gives reference so that we can test convergence of
% numerics, particularly in high-contrast case.)
% source: from H.Ammari, M.Putinar, 2016

% n: number of expansion terms. 
% rad: radius of circle
% epsilon1: inside permittivity.
% delta: separating distance between two circles. 

alpha = sqrt(delta*(rad + delta/4 ));

xi0 = asinh(alpha/rad);

lambda = (epsilon1 +1)/(2*(epsilon1 - 1));


num_t = 1:n;

vect = sum((num_t.*exp(-2.*num_t*xi0))./(lambda - 0.5*exp(-2.*num_t*xi0))); 

polariz = (8*pi*alpha^2)*vect;

end