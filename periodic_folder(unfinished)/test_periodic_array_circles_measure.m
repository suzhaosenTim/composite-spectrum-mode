clear




cparams_i = [];

cparams_i.eps = 1e-6;


rad = 0.4;
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
    weights = chnkrs(k).wts(:).';
    vals = vals + skern(chnkrs(k), chnkr).*weights;
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

[rho, lambda] = eig(lhs);

spec = diag(lambda);
iid = (spec > 1e-6);
spec_E = sort(spec(iid));


nx1 = chnkr.n(1,:).';

coef = rho\nx1;                         % Coefficient of e_1\cdot\nu in the eigenbasis





xi = -X*rho;
dens = rho + H*(R.'*xi);            % this step is somewhat redudant for 
                                    % eigenvectors we care about as they
                                    % are always orthogonal to R. Keep it
                                    % for consistency.

             
w = chnkr.wts(:);
xy = reshape(chnkr.r, 2, []);

measure_wts = (xy.* w.')* dens;        % wts formula integral of (\phi_n (e_1\cdot x) ds, \phi_n (e_2\cdot x) ds)
lambda = diag(lambda);

mom_0 = sum(measure_wts(1,:).'.*coef);
mom_1 = sum(lambda.*measure_wts(1,:).'.*coef);

epsilon1 = 10^8; epsilon2 = 1;

s = 1/(1-epsilon1/epsilon2);

vals_extr = epsilon2 - epsilon2*sum(1./(s-lambda).*measure_wts(1,:).'.*coef);


ref_sol = 3.10469989397989;

true_vals_1 = (pi*rad^2*(1-pi*rad^2))/2;            % first moment analytical formula p_1p_2/2


norm(mom_0 - pi*rad^2)
norm(mom_1 - true_vals_1)
norm(vals_extr - ref_sol)


wts_normal = (measure_wts(1,:).'.*coef);
[lambda_data, id] = sort(lambda);

CDF_mu = cumsum(wts_normal(id))/ sum(wts_normal);



figure(3)
plot(lambda_data(2:end),CDF_mu(2:end), 'LineWidth', 1.5)
grid on
xlabel('s')
ylabel('P((0,s])')

title('Cumulative Dstribution Function ')   

figure(4)
plot(1:length(spec_E), sort(spec_E), '.','MarkerSize', 18)
xlabel('index', 'FontSize',12)
ylabel('eigenvalues', 'FontSize',12)
axis square


