clear

cparams = [];
cparams.ta = 0;
cparams.tb = 1;
cparams.ifclosed = false;

e_1 = [1;0];
e_2 = [0;1];

chnkr1 = chunkerfunc(@(t) rotateline(t, 0, -0.5), cparams);
chnkr1 = chnkr1.reverse();

chnkr2 = chunkerfunc(@(t) rotateline(t, pi/2, -0.5), cparams);



M = 1:170;                          % controls the number of pxypoints;

max_err = zeros(1, length(M));
solv_err = zeros(1, length(M));

for k = 1:numel(M)

    j = 0:(M(k)-1);                 % setting parameters for sources

    R_p = 1.5;                          % setting pxypoints radius

    src_pts = [];

    src_pts.r = R_p.*[cos(2*pi*j./(M(k))); sin(2*pi*j./(M(k)))];



    fkern = @(s,t) chnk.lap2d.kern(s,t, 's');            % single layer kernel G(x,y)

    skern = @(s,t) chnk.lap2d.kern(s,t, 'sprime');        % \partial_nu G(x,y)

    Ded = chnkr1; Upp = chnkr1.move([], e_2); 

    Lef = chnkr2;  Rig= chnkr2.move([], e_1);


    lhs = generatePeriodiz(chnkr1,chnkr2, M(k), R_p);

    % analytical solution test with pt-src 1.5(cos(0.7); sin(0.7)) (same
    % with the paper)

    targs = [];
    targs.r = 1.5.*[cos(0.7); sin(0.7)];

    rhs1 =  fkern(targs, Rig) - fkern(targs, Lef);
    rhs2 =  skern(targs, Rig) - skern(targs, Lef);      % flip the sign due to the negative sign convention in chunkie
    rhs3 =  fkern(targs, Upp) - fkern(targs, Ded);
    rhs4 =  skern(targs, Upp) - skern(targs, Ded);     % flip the sign due to the negative sign convention in chunkie

    rhs = [rhs1; rhs2; rhs3; rhs4];

    % do a least-square solve via linsolve
    
    opts1 = [];
    opts1.RECT = true;

    coefs = linsolve(lhs, rhs, opts1); 

    solv_err(k) = norm(lhs*coefs - rhs);

    [xx, yy] = meshgrid(linspace(-.1,.1,10));
    eval_targs = []; eval_targs.r = [xx(:).'; yy(:).'];

    eval_mat = fkern(src_pts, eval_targs); 

    sol = eval_mat*coefs; 
   
    true_sol = fkern(targs, eval_targs);

    err = sol - true_sol;          % expect a constant generate by the null-space; 

    err = err - err(1);

    max_err(k) = max(err);
end

semilogy(M, max_err, 'LineStyle','-');
hold on 
loglog(M, 1./(M.^(7)), 'LineStyle','--')

figure(2); clf

loglog(M, solv_err)








