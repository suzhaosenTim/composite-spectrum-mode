clear

% the style of generating geometry is from the demo_mult_connect.m

iseed = 8675309;
rng(iseed,'twister');


cparams = [];
cparams.eps = 1.0e-6;
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



fkern = @(s,t) chnk.lap2d.kern(s,t,'sprime');

opts = []; 
opts.sing = 'log';
opts.adaptive_correction= true; % use adaptive quadrature to handle near touching panels.

start = tic;
sysmat = chunkermat(chnkrs, fkern, opts);
t1 = toc(start);
fprintf('%5.2e s : time to discretize the integral operator\n',t1)

epsi = 1e-3;
s = linspace(0,1, 200) ;

polariz_x = zeros(1, length(s));
polariz_y = zeros(1, length(s));
x1 = chnkrs.r(1,:);
x2 = chnkrs.r(2,:);
for i = 1:length(s)
    s1 = s(i) + epsi*1i;
    lhs = (s1- 0.5).*eye(chnkrs.npt) - sysmat;
    rhs = s1*chnkrs.n(1,:).';
    rhs2 = s1*chnkrs.n(2,:).';
    start1 = tic;
    rho = lhs\rhs;
    rho2 = lhs\rhs2;
    t2 = toc(start1);
    fprintf('%5.2e s : time for dense inversion\n',t2)
    
    
    polariz_x(i) = -(chnkrs.wts(:).')*(rho.*x1.')/s1;
    polariz_y(i) = -(chnkrs.wts(:).')*(rho2.*x2.')/s1;
end

fig = figure;
clf

t = tiledlayout(fig,1,2, ...
    'TileSpacing','compact', ...
    'Padding','loose');


nexttile
plot(real(s), real(polariz_x), 'LineWidth', 2)
hold on 
plot(real(s), imag(polariz_x), 'LineWidth',2)
xlabel('s-plane', 'FontSize',12)
ylabel('Polarizability', 'FontSize',12)
legend({'Real-part', 'Imaginary part'}, 'FontSize',12)
axis square
pbaspect([1 1 1]) 



nexttile
plot(real(s), real(polariz_y), 'LineWidth', 2)
hold on 
plot(real(s), imag(polariz_y), 'LineWidth',2)
xlabel('s-plane', 'FontSize',12)
ylabel('Polarizability', 'FontSize',12)
legend({'Real-part', 'Imaginary part'}, 'FontSize',12)
axis square
pbaspect([1 1 1]) 

set(fig, 'Units','pixels', 'Position',[50 50 1500 700]);