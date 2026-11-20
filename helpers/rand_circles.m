function [chnkr_int,ctrs, chnkrs] = rand_circles(L, rad, ninc)
% method of generating randomly placed circles within a fixed ''box''
% (code style followed from demo_mult_connect.m)

% L: size of the box we put the circles in
% rad: radius of the circle
% ninc: # of inclusions.



cparams = [];
cparams.eps = 1.0e-6;


chnkr_int = [];
ctrs = zeros(2,0);    
ntry = 1000;
for i = 1:ninc
    placed = false;    
    chnkr_i = chunkerfunc(@(t) rad .* circle(t), cparams);
    for j = 1:ntry
        tmp = (L - 2*rad)*rand(2,1);
        if isempty(ctrs)
            dist_ok = true;
        else
             dists = vecnorm(tmp - ctrs);
             dist_ok = all(dists > 2*rad); 
        end

        if dist_ok
           chnkr_i = chnkr_i.move([],tmp );
           chnkr_int = [chnkr_int, chnkr_i];
           ctrs(:,end+1) = tmp;
           placed = true;
           break
        end
    end

    if ~placed
        error('Could not place circle without touching, enlarge the wall L');
    end

end

chnkrs = merge(chnkr_int);
end