clear


iseed = 8675309;
rng(iseed,'twister');


L = 1.5;
rad = 0.1;
ninc = 16;

[~, ~, chnkrs] = rand_circles(L, rad, ninc);

plot(chnkrs);
hold on 
axis equal




fkern = @(s,t) chnk.lap2d.kern(s,t,'sprime');

opts = []; 
opts.sing = 'log';
opts.adaptive_correction= true; % use adaptive quadrature to handle near touching panels.

start = tic;
sysmat = chunkermat(chnkrs, fkern, opts);
t1 = toc(start);
fprintf('%5.2e s : time to discretize the integral operator\n',t1)


lhs = 0.5.*eye(chnkrs.npt) + sysmat  ;     % left hand side

[rho, lambda_max] = eigs(lhs, 1, 'largestreal');

spec = eig(lhs);
iid = (spec > 1e-6);
spec_E = spec(iid);







L = max(vecnorm(chnkrs.r(:,:)));

x1 = linspace(-0.2, 1.4 ,500);               % generate some targets to evaluate the eigenfunctions
[xx,yy] = meshgrid(x1,x1);
targets = [xx(:).'; yy(:).'];
                   

in = chunkerinterior(chnkrs, targets);              % identify interior and exterior (there is a hacky way to do this via FMM)
out = ~in;

eval_kern = @(s,t) chnk.lap2d.kern(s,t,'sgrad');  % gradient of the single layer

start2 = tic;
gamma_eig_in= chunkerkerneval(chnkrs, eval_kern, rho, targets(:,in));
gamma_eig_out = chunkerkerneval(chnkrs, eval_kern, rho, targets(:,out));
t3 = toc(start2);
fprintf('%5.2e s : time to plot \n',t3)

gamma_eigvec = nan(2,size(xx(:),1));

gamma_eigvec(:,in(:)) = reshape(gamma_eig_in, 2, []);

gamma_eigvec(:,out(:)) = reshape(gamma_eig_out, 2,[]);

strE = abs(gamma_eigvec(:, in(:)));

strE_out = abs(gamma_eigvec(:, out(:)));

str_plot = nan(size(xx));

str_plot(in) = sqrt(strE(1,:).^2 + strE(2,:).^2); 
str_plot(out) = sqrt(strE_out(1,:).^2 + strE_out(2,:).^2); 



fig = figure;
clf

t = tiledlayout(fig,1,3, ...
    'TileSpacing','compact', ...
    'Padding','loose');


nexttile
h = pcolor(xx, yy, (str_plot));
set(h,'EdgeColor','none');
hold on
plot(chnkrs,'k')
xlabel('x-axis', 'FontSize',12)
ylabel('y-axis', 'FontSize',12)
%title('Eigenfunction (Linear scaled absolute value) ', 'FontSize',12)
axis equal
colorbar
pbaspect([1 1 1]) 

nexttile
h = pcolor(xx, yy, log10(str_plot));
set(h,'EdgeColor','none');
hold on
plot(chnkrs,'k')
xlabel('x-axis', 'FontSize',12)
ylabel('y-axis', 'FontSize',12)
%title('Eigenfunction (Log scaled absolute value)', 'FontSize',12)
axis equal
colorbar
pbaspect([1 1 1]) 

nexttile
plot(1:length(spec_E), sort(spec_E), '.','MarkerSize', 18)
xlabel('index', 'FontSize',12)
ylabel('eigenvalues', 'FontSize',12)
axis square
pbaspect([1 1 1]) 

set(fig, 'Units','pixels', 'Position',[50 50 1500 700]);



% 
exportgraphics(gcf,'null_vector_unordered.pdf','ContentType','vector',...
   'Resolution',600)
