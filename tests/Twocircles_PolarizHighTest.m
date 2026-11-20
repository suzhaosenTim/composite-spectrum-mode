Twocircles_polarizTest0();


function Twocircles_polarizTest0()


cparams = [];

cparams.eps = 1e-6;
cparams.nover = 0;


rad = 0.1;
delta = 0.256;                        % separation
ctr = 0;

% two circle case
start = tic; 
chnkrorg = chunkerfunc(@(t) rad.*circle(t),cparams);
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)
    
ctr1 = [ctr + delta/2; 0];
ctr2 = [-ctr - delta/2; 0];

chnkr = chnkrorg.move([], ctr1 );

chnkr1 = chnkrorg.move([], ctr2);

chnkr_int = [chnkr, chnkr1];

chnkrs = merge(chnkr_int);


fkern = @(s,t) chnk.lap2d.kern(s,t,'sprime');

opts = []; 
opts.sing = 'log';

start = tic;
sysmat = chunkermat(chnkrs, fkern, opts);
t1 = toc(start);
fprintf('%5.2e s : time to discretize the integral operator\n',t1)

epsilon1 = 10^(15);   epsilon2 = 1;                 % inside permittivity 10^15

s = 1/(1-epsilon1/epsilon2);


stb = stb_multiple(chnkr_int);    
    
lhs = (0.5-s).*eye(chnkrs.npt) + sysmat + stb ;

rhs = -chnkrs.n(1,:).';
rho = lhs\rhs;

x1 = chnkrs.r(1,:);
wts = chnkrs.wts;

integrand = rho.*x1.';

polariz_num = -(wts(:).')*integrand;

n = 400;
delta1 = delta -2*rad ;

polariz = asymptot_high_contrast(n,rad, epsilon1, delta1);      % asymptotic expansion of High contrast polarizability

rel_err = norm(polariz-polariz_num)/(norm(polariz));
abs_err = abs(polariz - polariz_num);

assert(rel_err < 1e-8)
assert(abs_err <1e-8)

end