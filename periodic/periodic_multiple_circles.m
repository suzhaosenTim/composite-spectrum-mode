clear

iseed = 8675309;
rng(iseed,'twister');


L = 1.0;
rad = 0.07;
ninc = 16;

[~, ~, chnkrorg] = rand_circles(L, rad, ninc);

chnkrs = chnkrorg.move([], [-0.4;-0.4]);


cparams = [];
cparams.ta = 0;
cparams.tb = 1.0;
cparams.ifclosed = false;

start = tic; 
chnkr1 = chunkerfunc(@(t) rotateline(t, 0, -0.5), cparams); % cell chnkrs
chnkr1 = chnkr1.reverse();                                  
% make sure the normal direction aligned with the convention.
chnkr2 = chunkerfunc(@(t) rotateline(t, pi/2, -0.5), cparams);
t2 = toc(start);
fprintf('%5.2e s : time to build cells\n',t2)

% kernels
fkern = @(s,t) chnk.lap2d.kern(s,t, 's');
skern = @(s,t) chnk.lap2d.kern(s,t, 'sprime');

opts = [];
opts.sing = 'log';
opts.adaptive_correction = true;

start2 = tic;
mat1 = chunkermat(chnkrs, skern, opts);
t3 = toc(start2);
fprintf('%5.2e s : time to discretize self-interaction matrix\n',t3)


Rp = 3;
M = 100; 


[x, y] = meshgrid(-1:1); 
tral = [x(:).'; y(:).']; 
tral = [tral(:,1:4), tral(:,6:end)];


chnkrs_img =  sumSrccopies(chnkrs, tral); 


vals = 0;
for k = 1:numel(chnkrs_img)
    vals = vals + chunkerkernevalmat(chnkrs, skern, chnkrs_img(k));
end

A_block = 0.5.*eye(chnkrs.npt) + mat1 + vals;



j = 0:(M-1);
src_pts = [];
src_pts.r = Rp.*[cos(2*pi*j./(M)); sin(2*pi*j./(M))];              % pxypts used to form block B

Bmat = skern(src_pts, chnkrs);                                       % (original B matrix)

N = chnkrs.npt;

R = 1/M.*ones(M, 1); H = ones(N, 1);

B = Bmat + A_block*H*R.';                                       % new B_block for Schur complement


ctr = [0;0];

start5 = tic;
[objcchnkr,Cmat] = formCblockmat(chnkrs,ctr, fkern, skern, chnkr1, chnkr2);

Qmat = generatePeriodiz(chnkr1, chnkr2, M, Rp); 
t5 = toc(start5);
fprintf('%5.2e s : time to form wall matrix C and Q\n',t5)


Q = Qmat + Cmat*H*R.';

X = Q\Cmat;

lhs = A_block - B*X;

[rho, lambda] = eig(lhs);

s2 = diag(lambda);
iid = (s2 > 1e-6);      
spec_E = sort(s2(iid));


nx1 = chnkrs.n(1,:).';

coef = rho\nx1;                         % Coefficient of e_1\cdot\nu in 
                                        % the eigenbasis

xi = -X*rho;
dens = rho + H*(R.'*xi);            % this step is somewhat redudant for 
                                    % eigenvectors we care about as they
                                    % are always orthogonal to R. Keep it
                                    % for consistency.

             
w = chnkrs.wts(:);
xy = reshape(chnkrs.r, 2, []);

measure_wts = (xy.* w.')* dens;        % wts formula integral of (\phi_n (e_1\cdot x) ds, \phi_n (e_2\cdot x) ds)


mom_0 = sum(measure_wts(1,:).'.*coef);
true_ans = pi*rad^2*ninc/(L^2);

abs_err = norm(mom_0 - true_ans);
rel_err = norm(mom_0 - true_ans)/norm(true_ans);


assert(abs_err < 1e-8)
assert(rel_err < 1e-8)