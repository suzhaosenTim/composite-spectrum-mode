function Qmat = generatePeriodiz(chnkr1, chnkr2, M, Rp)
% A method to generate the periodized matrix Q\xi = sum_{j = 1}^M xi_j
% G(x,y_j), where G is the laplace 2d kernel.
% 
% See A.Barnett 2018 CPAM paper
%
% Input two chunker object (for unit cell. Typically for squared cells)
% M: # of pxypoints
% Rp: radius of the pxypoints

j = 0:(M-1); 
src_pts = [];
src_pts.r = Rp.*[cos(2*pi*j./(M)); sin(2*pi*j./(M))];       % proxy points


fkern = @(s,t) chnk.lap2d.kern(s,t, 's');            % single layer kernel G(x,y)

skern = @(s,t) chnk.lap2d.kern(s,t, 'sprime');        % \partial_nu_x G(x,y)

e_2 = [0;1];
e_1 = [1;0];

Ded = chnkr1; 
Upp = chnkr1.move([], e_2); 

Lef = chnkr2; 
Rig = chnkr2.move([], e_1);


Q1mat = fkern(src_pts, Rig) - fkern(src_pts, Lef);
Q2mat = skern(src_pts, Rig) - skern(src_pts, Lef);   
Q3mat = fkern(src_pts, Upp) - fkern(src_pts, Ded);
Q4mat = skern(src_pts, Upp) - skern(src_pts, Ded);  

Qmat = [Q1mat; Q2mat;Q3mat; Q4mat];

end

