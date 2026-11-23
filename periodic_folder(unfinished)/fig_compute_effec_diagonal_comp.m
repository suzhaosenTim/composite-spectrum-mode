clear

iseed = 8675309;
rng(iseed,'twister');



rad = 0.1;
delta = 0.256;                        % separation
ctr = 0;

start = tic; 
chnkrorg = chunkerfunc(@(t) rad.*circle(t));
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)
    
ctr1 = [ctr + delta/2; 0];
ctr2 = [-ctr - delta/2; 0];

chnkr = chnkrorg.move([], ctr1 );
chnkr1 = chnkrorg.move([], ctr2);

chnkr_int = [chnkr, chnkr1];

chnkrs = merge(chnkr_int);



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
iid = (s2>1e-6);
spec_E = s2(iid);



nx1 = chnkrs.n(1,:).';
nx2 = chnkrs.n(2,:).';

coef = rho\nx1;                         % Coefficient of e_1\cdot\nu in 
                                        % the eigenbasis
coef_2 = rho\nx2;
                                      

xi = -X*rho;
dens = rho + H*(R.'*xi);            % this step is somewhat redudant for 
                                    % eigenvectors we care about as they
                                    % are always orthogonal to R. Keep it
                                    % for consistency.

             
w = chnkrs.wts(:);
xy = reshape(chnkrs.r, 2, []);

measure_wts = (xy.* w.')* dens;        % wts formula integral of (\phi_n (e_1\cdot x) ds, \phi_n (e_2\cdot x) ds)

mom_0 = sum(measure_wts(1,:).'.*coef);

true_ans = pi*rad^2*2;
abs(mom_0 - true_ans)


s = linspace(0,1 , 5000) + 1e-3i;

vals_extr = 1 - sum(1./(s-s2).*measure_wts(1,:).'.*coef);
vals_extr_2 = 1 - sum(1./(s-s2).*measure_wts(2,:).'.*coef_2);

fig = figure;
clf

t = tiledlayout(fig,1,3, ...
    'TileSpacing','compact', ...
    'Padding','loose');
nexttile
plot(real(s), real(vals_extr), 'LineWidth',2);
hold on 
plot(real(s), imag(vals_extr),'LineWidth',2);
xlabel('s-plane', 'FontSize',12)
ylabel('Effective permittivity (11)', 'FontSize',12)
legend({'Real-part', 'Imaginary part'}, 'FontSize',12)
axis square
pbaspect([1 1 1]) 

nexttile
plot(real(s), real(vals_extr_2), 'LineWidth',2);
hold on 
plot(real(s), imag(vals_extr_2),'LineWidth',2);
xlabel('s-plane', 'FontSize',12)
ylabel('Effective permittivity (22)', 'FontSize',12)
legend({'Real-part', 'Imaginary part'}, 'FontSize',12, 'Location','northwest')
axis square
pbaspect([1 1 1]) 

nexttile
plot(1:length(spec_E), sort(spec_E),'.', 'MarkerSize',18);
xlabel('index', 'FontSize',12)
ylabel('eigenvalues', 'FontSize',12)
axis square
pbaspect([1 1 1]) 

set(fig, 'Units','pixels', 'Position',[50 50 1500 700]);

exportgraphics(gcf,'Effective permittivity as a function of s.pdf','ContentType','vector')