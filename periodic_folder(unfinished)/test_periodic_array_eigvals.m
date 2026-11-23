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
% make sure the normal direction aligned with the convention.
chnkr2 = chunkerfunc(@(t) rotateline(t, pi/2, -0.5), cparams);
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

[rho, lambda] = eigs(lhs,10, 'smallestreal');



rho = rho(:,2);
xi = -X*rho;        % pxy coefficient

dens = rho + H*R.'*xi;             % this step is somewhat redudant for eigenvectors 
                                    % as the ones we care about
                                    % are always orthogonal to R. Keep it
                                    % for consistency.

 
plot(objcchnkr); hold on;   scatter(src_pts.r(1,:), src_pts.r(2,:),'*');
axis equal;

x1 = linspace(-1.4, 1.4  ,500);
[xx,yy] = meshgrid(x1,x1);
targets = [xx(:).'; yy(:).'];




tic
in = chunkerinterior(objcchnkr, targets) | chunkerinterior(chnkr, targets);
out = ~in;
toc

% grad of single layer part.
gkern = @(s,t) chnk.lap2d.kern(s,t , 'sgrad');

start6 = tic;
Dsol = chunkerkerneval(chnkr, gkern, dens, targets(:,out)); 

for k = 1:length(chnkrs)
    Dsol = Dsol + chunkerkerneval(chnkrs(k), gkern, dens, targets(:,out));
end

Dsol_i = chunkerkerneval(chnkr, gkern, dens, targets(:,in)); 

for k = 1:length(chnkrs)
    Dsol_i = Dsol_i + chunkerkerneval(chnkrs(k), gkern, dens, targets(:,in));
end
t6 = toc(start6);
fprintf('%5.2e s : time to evaluate eigenfunctions\n',t6)

Dsol = reshape(Dsol, 2, []);
Dsol_i = reshape(Dsol_i, 2,[]);



% pxy_points part

targs = []; 
targs.r = targets(:,out);
targs_i = []; 
targs_i.r = targets(:,in);

eval_mat = gkern(src_pts, targs); 
periodicsol_out = eval_mat*xi;

eval_mat_i = gkern(src_pts, targs_i);
periodicsol_in = eval_mat_i*xi;


% combine them 
sol = Dsol + reshape(periodicsol_out, 2,[]);
sol_i = Dsol_i + reshape(periodicsol_in, 2, []);


uugrad = nan(2,size(xx(:),1));
uugrad(:, in(:)) = sol_i;
uugrad(:, out(:)) = sol;

strE = reshape(vecnorm(uugrad, 2, 1), size(xx));

figure(1)
clf
h = pcolor(xx, yy, strE);
set(h,'EdgeColor','None'); hold on;
plot(chnkrs,'k','LineWidth',1);
hold on 
plot(chnkr, 'k','LineWidth',1);
xlabel('x-axis', 'FontSize',12)
ylabel('y-axis', 'FontSize',12)
title('Eigfunc ', 'FontSize',12)
colorbar
axis equal

figure(2)
clf
h = pcolor(xx, yy, log10(strE));
set(h,'EdgeColor','None'); hold on;
plot(chnkrs,'k','LineWidth',1);
hold on 
plot(chnkr, 'k','LineWidth',1);
xlabel('x-axis', 'FontSize',12)
ylabel('y-axis', 'FontSize',12)
title('Eigfunc ', 'FontSize',12)
colorbar
axis equal
