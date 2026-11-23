clear



cparams = [];

cparams.eps = 1e-6;
cparams.nover = 0;
%cparams.maxchunklen = 0.5;

rad = 0.1;
delta = 0.256;
ctr = 0;


nch = 1:1:8;

rel_err = zeros(1, length(nch));
npts = zeros(1, length(nch));
cond_number = zeros(1, length(nch));

for ii = 1: length(rel_err)
    start = tic; 
    chnkrorg = chunkerfuncuni(@(t) rad.*circle(t),nch(ii), cparams); 
    t1 = toc(start);
    fprintf('%5.2e s : time to build geo\n',t1)

    ctr1 = [ctr + delta/2; 0];

    ctr2 = [-ctr - delta/2; 0];

    chnkr = chnkrorg.move([], ctr1 );

    chnkr1 = chnkrorg.move([], ctr2);

    ctrs = [ctr1, ctr2];
 


    chnkr_int = [chnkr, chnkr1];

    chnkrs = merge(chnkr_int);
    npts(ii)= chnkrs.npt;


    fkern = @(s,t) chnk.lap2d.kern(s,t,'sprime');


    opts = []; 
    opts.sing = 'log';

    start = tic;
    sysmat = chunkermat(chnkrs, fkern, opts);
    t1 = toc(start);
    fprintf('%5.2e s : time to discretize the integral operator\n',t1)

    epsilon1 = 10^(15);   epsilon2 = 1;

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
    delta1 = delta - 2*rad;

    polariz = asymptot_high_contrast(n,rad, epsilon1, delta1);

    rel_err(ii) = norm(polariz-polariz_num)/(norm(polariz));
end


fig = figure;
clf

t = tiledlayout(fig,1,2, ...
    'TileSpacing','compact', ...
    'Padding','loose');

nexttile
loglog(npts(1:8), rel_err(1:8), '.-', 'LineWidth',2,'MarkerSize',12)
hold on 
loglog(npts(1:8), (10^22).*npts(1:8).^(-14), '--k', 'LineWidth',1.5)
grid on
legend({'polarizability relative error', 'n^{-14}'}, 'FontSize',12)
xlabel('Discretization nodes (Total)', 'FontSize',12)
ylabel('Relative error','FontSize',12)


axis square



nexttile
plot(chnkrs,'-k', 'LineWidth',2)
hold on 
plot(chnkrs, 'xk', 'LineWidth',2)
axis equal
pbaspect([1 1 1]) 

set(fig, 'Units','pixels', 'Position',[50 50 1000 700]);


exportgraphics(gcf,'convergence_figure_polarizability.pdf','ContentType','vector')