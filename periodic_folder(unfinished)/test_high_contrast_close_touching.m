clear




cparams_i = [];

cparams_i.eps = 1e-6;


rad = 0.5 - 1e-6;
start = tic; 
chnkr = chunkerfunc(@(t) rad.*circle(t),cparams_i);
%chnkr = refine(chnkr, struct('nover', 1));

t1 = toc(start);
fprintf('%5.2e s : time to build geo\n',t1)


    

cparams = [];
cparams.ta = 0;
cparams.tb = 1;
cparams.ifclosed = false;

start = tic; 
chnkr1 = chunkerfunc(@(t) rotateline(t, 0, -0.5), cparams); % cell chnkrs
chnkr1 = chnkr1.reverse();
%chnkr1 = refine(chnkr1, struct('nover', 1));
% make sure the normal direction aligned with the convention.
chnkr2 = chunkerfunc(@(t) rotateline(t, pi/2, -0.5), cparams);
%chnkr2 = refine(chnkr2, struct('nover', 1));
t2 = toc(start);
fprintf('%5.2e s : time to build cells\n',t2)


% kernels
fkern = @(s,t) chnk.lap2d.kern(s,t, 's');
skern = @(s,t) chnk.lap2d.kern(s,t, 'sprime');

opts = [];
opts.sing = 'log';
opts.adaptive_correction = true;

start2 = tic;
mat1 = chunkermat(chnkr, skern, opts);
t3 = toc(start2);
fprintf('%5.2e s : time to discretize self-interaction matrix\n',t3)


Rp = 3;
M = 100; 




[x, y] = meshgrid(-1:1); 
tral = [x(:).'; y(:).']; 
tral = [tral(:,1:4), tral(:,6:end)];


chnkrs =  sumSrccopies(chnkr, tral); 


vals = 0;
for k = 1:numel(chnkrs)
    vals = vals + chunkerkernevalmat(chnkr, skern, chnkrs(k));
end

A_block = 0.5.*eye(chnkr.npt) + mat1 + vals;



j = 0:(M-1);
src_pts = [];
src_pts.r = Rp.*[cos(2*pi*j./(M)); sin(2*pi*j./(M))];              % pxypts used to form block B

Bmat = skern(src_pts, chnkr);                                       % (original B matrix)

N = chnkr.npt;

R = 1/M.*ones(M, 1); H = ones(N, 1);

B = Bmat + A_block*H*R.';                                       % new B_block for Schur complement


ctr = [0;0];

start5 = tic;
[objcchnkr,Cmat] = formCblockmat(chnkr,ctr, fkern, skern, chnkr1, chnkr2);

Qmat = generatePeriodiz(chnkr1, chnkr2, M, Rp); 
t5 = toc(start5);
fprintf('%5.2e s : time to form wall matrix C and Q\n',t5)


Q = Qmat + Cmat*H*R.';

X = Q\Cmat;

lhs = A_block - B*X;

[~, lambda] = eig(lhs);

s2 = diag(lambda);
max(s2)

scatter(1:length(s2), s2, "*")

%%
epsilon1 = 10^3; epsilon2 = 1;                             % permittivity

s = 1/(1- epsilon1/epsilon2);

Aself = (0.5-s).*eye(chnkr.npt) + mat1 + onesmat(chnkr); 
lhs = [Aself Bmat; 
       Cmat Qmat];                                          % form the rectangular lhs


rhs1 = zeros(chnkr.npt, 1);                     
rhs2 = [ones(chnkr1.npt, 1);                                % m = 16 (we use # 16 nodes)
        zeros(chnkr1.npt, 1);
        zeros(chnkr1.npt, 1);
        zeros(chnkr1.npt, 1)];


rhs = [rhs1; rhs2];



lso.RECT = true;

tic
dens = linsolve(lhs, rhs, lso);
toc

rho = dens(1:chnkr.npt);                        % integral density
xi = dens((chnkr.npt+ 1):end);                  % pxy density

 
x1 = chnkr.r(1,:).';
epsilon_eff_1 = epsilon2 - epsilon2*chunkerintegral(chnkr, rho.*x1)