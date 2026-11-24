clear




cparams = [];
cparams.ta = 0;
cparams.tb = 1;
cparams.ifclosed = false;

e_1 = [1;0];
e_2 = [0;1];

chnkr1 = chunkerfunc(@(t) rotateline(t, 0, -0.5), cparams);
chnkr1 = chnkr1.reverse();
chnkr2 = chnkr1.move([], e_2);

chnkr3 = chunkerfunc(@(t) rotateline(t, pi/2, -0.5), cparams);


chnkr4 = chnkr3.move([], e_1);

chnkr_int = [chnkr1, chnkr2, chnkr3, chnkr4];

chnkrs = merge(chnkr_int);
plot(chnkrs)
quiver(chnkrs)
axis equal
% hold on 
% plot(chnkr4)
% quiver(chnkr4)