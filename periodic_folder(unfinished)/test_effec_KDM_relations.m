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
cparams.tb = 1;
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


j = 0:(M-1);
src_pts = [];
src_pts.r = Rp.*[cos(2*pi*j./(M)); sin(2*pi*j./(M))];              % pxypts used to form block B

Bmat = skern(src_pts, chnkrs);                                       % (original B matrix)

N = chnkrs.npt;

R = 1/M.*ones(M, 1); H = ones(N, 1);

ctr = [0;0];

start5 = tic;
[objcchnkr,Cmat] = formCblockmat(chnkrs,ctr, fkern, skern, chnkr1, chnkr2);

Qmat = generatePeriodiz(chnkr1, chnkr2, M, Rp); 
t5 = toc(start5);
fprintf('%5.2e s : time to form wall matrix C and Q\n',t5)


Q = Qmat + Cmat*H*R.';

X = Q\Cmat;

h = 0.8; s = 1/(1-h);

h2 = 1/h; s2 = 1/(1-h2);

A_block_1 = (s-0.5).*eye(chnkrs.npt) + mat1 + vals;
A_block_2 = (s2-0.5).*eye(chnkrs.npt) + mat1 + vals;

lhs = A_block_1 - Bmat*X;
lhs2 = A_block_2 - Bmat*X;

                     
rhs = s.*[ones(chnkr1.npt, 1);                                
        zeros(chnkr1.npt, 1);
        zeros(chnkr1.npt, 1);
        zeros(chnkr1.npt, 1)];

rhs2 = s.*[zeros(chnkr1.npt, 1);                                
        zeros(chnkr1.npt, 1);
        ones(chnkr1.npt, 1);
        zeros(chnkr1.npt, 1)];

rhs3 = s2.*[ones(chnkr1.npt, 1);                                
        zeros(chnkr1.npt, 1);
        zeros(chnkr1.npt, 1);
        zeros(chnkr1.npt, 1)];

rhs4 =  s2.*[zeros(chnkr1.npt, 1);                                
        zeros(chnkr1.npt, 1);
        ones(chnkr1.npt, 1);
        zeros(chnkr1.npt, 1)];

rhs_h_1 = -Bmat*(Q\rhs);
rhs_h_2 = -Bmat*(Q\rhs2);

rhs_h2_1 = -Bmat*(Q\rhs3);
rhs_h2_2 = -Bmat*(Q\rhs4);

rho1 = lhs\rhs_h_1;
rho2 = lhs\rhs_h_2;

rho3 = lhs2\rhs_h2_1;
rho4 = lhs2\rhs_h2_2;

x1 = chnkrs.r(1,:);
x2 = chnkrs.r(2,:);

w = chnkrs.wts(:);

effec12 = -(x2.*w.')*rho1/s;
effec21 = -(x1.*w.')*rho2/s;

effec_2_12 = - (x2.*w.')*rho3/s2;
effec_2_21 = - (x1.*w.')*rho4/s2;

effec_2_11 = 1 - (x1.*w.')*rho3/s2;
effec_2_22 = 1 - (x2.*w.')*rho4/s2;



% check Keller-Dykhne-Menedelson relations


effec_2 = [effec_2_11, effec_2_12;
           effec_2_21, effec_2_22];

effec_inv = effec_2\eye(2);


norm(effec_inv(1,2) - (-effec12))
norm(-effec_inv(2,1 ) - effec21 )











