clear

% same test problem as the one in the paper (quick implementation by
% solving the ELS via least squares) 



chnkr = chunkerfunc(@(t) wormcurv(t, 0.7, 0.15));   % inclusion chunk.
chnkr = refine(chnkr, struct('nover', 1));




cparams = [];
cparams.ta = 0;
cparams.tb = 1;
cparams.ifclosed = false;


chnkr1 = chunkerfunc(@(t) rotateline(t, 0, -0.5), cparams); % cell chnkrs
chnkr1 = chnkr1.reverse();                                  % make sure the normal direction.

chnkr2 = chunkerfunc(@(t) rotateline(t, pi/2, -0.5), cparams);


fkern = @(s,t) chnk.lap2d.kern(s,t, 's');
skern = @(s,t) chnk.lap2d.kern(s,t, 'sprime');

opts = [];
opts.sing = 'log';

tic
mat1 = chunkermat(chnkr, skern, opts);
toc

Aself = -0.5.*eye(chnkr.npt) + mat1;                                 % self-interaction A

Rp = 1.5;
M = 70; 


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
Qmat = generatePeriodiz(chnkr1,chnkr2, M, Rp); 
toc

lhs = [Amat Bmat; 
       Cmat Qmat];                                          % form the rectangular lhs

[x1, y1] = meshgrid(-2:2); 
tral1 = [x1(:).'; y1(:).'];  
src_loc = [0; 0];                                         % 5x 5 source location

dkern = @(s,t) chnk.lap2d.kern(s,t, 'd');                   % D kernel
hkern = @(s,t) chnk.lap2d.kern(s, t, 'dprime');             % D' kernel 

src_images = [];

src_images.r = (src_loc + tral1);

src_images.n = [1;0];

rhsvec = hkern(src_images, chnkr);                          % check 

rhs1 = sum(rhsvec, 2);

rhs2 = formCblockrhs(src_loc, dkern, hkern, chnkr1,chnkr2, 2);   % build up the known discrepancies

rhs = [rhs1; rhs2];



lso.RECT = true;

tic
dens = linsolve(lhs, rhs, lso);
%dens = lhs\rhs;
toc



dens1 = dens(1:chnkr.npt); 
coefs = dens((chnkr.npt+1):end);           % extract densities of layer potentials and periodizing part




[xx, yy] = meshgrid(linspace(-.43,.43,10));

targets = [xx(:).'; yy(:).']; 


tic
in = chunkerinterior(objcchnkr, targets) | chunkerinterior(chnkr, targets);
out = ~in;
toc


targs = []; 
targs.r = targets(:, out);

eval_mat = fkern(src_pts, targs); 
periodicsol = eval_mat*coefs;




tic
Dsol = chunkerkerneval(chnkr, fkern, dens1, targets(:,out)); 
toc

for k = 1:length(chnkrs)
    Dsol = Dsol + chunkerkerneval(chnkrs(k), fkern, dens1, targets(:,out));
end



sol = Dsol + periodicsol; 

true_sol = dkern(src_images, targs); 

true_sol = sum(true_sol, 2); 

err = sol - true_sol; % should be off by a constant

err = abs(err - err(1));

max(err)                % max absolute error


