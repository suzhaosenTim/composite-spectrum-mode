clear




cparams_i = [];

cparams_i.eps = 1e-6;


rad = 0.4;
start = tic; 
chnkr = chunkerfunc(@(t) rad.*circle(t),cparams_i);

t1 = toc(start);
fprintf('%5.2e s : time to build geo\n',t1)


    

plot(chnkr)
axis equal


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

fkern = @(s,t) chnk.lap2d.kern(s,t, 's');
skern = @(s,t) chnk.lap2d.kern(s,t, 'sprime');

opts = [];
opts.sing = 'log';


tic
mat1 = chunkermat(chnkr, skern, opts);
toc

epsilon1 = 10^8; epsilon2 = 1;                             % permittivity

s = 1/(1- epsilon1/epsilon2);

Aself = (0.5-s).*eye(chnkr.npt) + mat1 + onesmat(chnkr);                                 % self-interaction A

Rp = 2;
M = 100; 


j = 0:(M-1);
src_pts = [];
src_pts.r = Rp.*[cos(2*pi*j./(M)); sin(2*pi*j./(M))];              % pxypts used to form block B

Bmat = skern(src_pts, chnkr); 

[x, y] = meshgrid(-1:1); 
tral = [x(:).'; y(:).']; 
tral = [tral(:,1:4), tral(:,6:end)];


chnkrs =  sumSrccopies(chnkr, tral); 


vals = 0;
for k = 1:numel(chnkrs)
    weights = chnkrs(k).wts(:).';
    vals = vals + skern(chnkrs(k), chnkr).*weights;
end


Amat = Aself + vals;                                                   % full block A

ctr = [0;0];
tic
[objcchnkr,Cmat] = formCblockmat(chnkr,ctr, fkern, skern, chnkr1, chnkr2);
toc

plot(objcchnkr); hold on;   scatter(src_pts.r(1,:), src_pts.r(2,:),'*');
axis equal;

tic
Qmat = generatePeriodiz(chnkr1, chnkr2, M, Rp); 
toc

lhs = [Amat Bmat; 
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
epsilon_eff_1 = epsilon2 - epsilon2*chunkerintegral(chnkr, rho.*x1);

ref_sol = 3.10469989397989;
norm(ref_sol - epsilon_eff_1)


L = max(abs(objcchnkr.r),[],"all");

x1 = linspace(-L, L  ,500);
[xx,yy] = meshgrid(x1,x1);
targets = [xx(:).'; yy(:).'];




tic
in = chunkerinterior(objcchnkr, targets) | chunkerinterior(chnkr, targets);
out = ~in;
toc


targs = []; 
targs.r = targets(:,out);

eval_mat = fkern(src_pts, targs); 
periodicsol_out = eval_mat*xi;

targs_i = []; 
targs_i.r = targets(:,in);
eval_mat_i = fkern(src_pts, targs_i); 
periodicsol_in = eval_mat_i*xi;

tic
Dsol = chunkerkerneval(chnkr, fkern, rho, targets(:,out)); 
toc

for k = 1:length(chnkrs)
    Dsol = Dsol + chunkerkerneval(chnkrs(k), fkern, rho, targets(:,out));
end


tic
Dsol_i = chunkerkerneval(chnkr, fkern, rho, targets(:,in)); 
toc

for k = 1:length(chnkrs)
    Dsol_i = Dsol_i + chunkerkerneval(chnkrs(k), fkern, rho, targets(:,in));
end

sol = Dsol + periodicsol_out;
sol_i = Dsol_i + periodicsol_in;

uupotential = nan(size(xx));

uupotential(in) = sol_i;
uupotential(out) = sol;

figure(3)
contour(xx, yy, abs(uupotential), 'LineWidth',1.5)
hold on 
plot(chnkrs,'-k')
title('Contour plot of potential', 'FontSize',12)
ylabel('potential','FontSize',12)
colorbar
axis equal


figure(4)
clf
h = pcolor(xx,yy,abs(uupotential));
set(h,'EdgeColor','None'); hold on;
plot(chnkrs,'k-','LineWidth',1);
xlabel('x-axis', 'FontSize',12)
ylabel('y-axis', 'FontSize',12)
title('|E| (log_{10} scale) ', 'FontSize',12)
axis equal
axis image
colorbar
