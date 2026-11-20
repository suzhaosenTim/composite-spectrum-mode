Twocircles_eigenvaluesTest0();


function Twocircles_eigenvaluesTest0()
% test eigenvalues of \Gamma \chi operator on a pair of circles. 
% analytic formula comes from E.Bonnetier and F.Triki's paper. 


cparams = [];

cparams.eps = 1e-6;
cparams.nover = 0;


rad = 0.1;
delta = 0.256;                        % separation
ctr = 0;

start = tic; 
chnkrorg = chunkerfunc(@(t) rad.*circle(t), cparams);
t1 = toc(start);
fprintf('%5.2e s : time to build geo\n',t1)
    
ctr1 = [ctr + delta/2; 0];
ctr2 = [-ctr - delta/2; 0];

chnkr = chnkrorg.move([], ctr1 );

chnkr1 = chnkrorg.move([], ctr2);

chnkr_int = [chnkr, chnkr1];

chnkrs = merge(chnkr_int);
chnkrs = refine(chnkrs, struct("nover", 1));        % refine 1-lvl     

fkern = @(s,t) chnk.lap2d.kern(s, t, 'sprime');

opts = [];
opts.sing = 'log';


start = tic;
sysmat = chunkermat(chnkrs, fkern, opts);
t1 = toc(start);
fprintf('%5.2e s : time to discretize the integral operator\n',t1)


lhs = 0.5.*eye(chnkrs.npt) + sysmat  ;     % left hand side

lambda = eig(lhs);

iid = lambda>1e-5;

lambda_corr = lambda(iid);                  % throw out the known 0 eigenvalue 

lambda_sing = lambda_corr(1:2:end);         % every eigenvalue has multiplicity 2

[n, ~] = size(lambda_sing);

% Compute two families max eigenvalues via formula in Bonnetier paper
% (slightly more general)

a = sqrt((ctr + (delta/2)-rad)*((ctr + (delta/2)-rad) + 2*rad));

rhodist = (ctr + (delta/2) -a)/rad ;


d2 = 1:n;

s_para_lower =  0.5 - rhodist.^(2.*d2)/2;
s_para_higher = 0.5 + rhodist.^(2.*d2)/2;

lambda_sing = sort(lambda_sing);

iid2 = (abs(lambda_sing - 0.5)>1e-5);               % throw out the trivial 1/2 accumulation point
lambda_sing = lambda_sing(iid2);

tru_eig = sort([s_para_lower, s_para_higher]);
iid3 = (abs(tru_eig - 0.5)>1e-5);
tru_eig = tru_eig(iid3);


abs_err = abs(lambda_sing - tru_eig.');
rel_err1 = norm(lambda_sing - tru_eig.')/norm(tru_eig.');

assert(max(abs_err)<1e-8)
assert(max(rel_err1)<1e-8)

end


