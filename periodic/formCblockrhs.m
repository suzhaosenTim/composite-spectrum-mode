function Cblockrhs = formCblockrhs(src_loc, kern1, kern2, chnkr1,chnkr2, n)
% A method to form the C Block rhs (only use for the analytic solution test);
% Input argument: chnkobj- a chunker or chunkgraph describing the inclusion
% ctr: center of the object
% kern1: first kernel- typically \partial_{n_y}G(x,y) 
% kern2: second kernel- typically \partial_n_x \partial_{n_y}G(x,y)
% (hypersingular kernel)

% Note: self-evaluation is handled separately. 

    [x, y] = meshgrid(-n:n); tral = [x(:).'; y(:).']; 

    i = length(x);                      % index to find the leftmost & rightmost copies
    
    paralef = [];
    pararig = [];
    paraded = [];
    paraupp = [];

    paralef.r = src_loc + tral(:,1:i);       % form the left-hand side i copies
    pararig.r = src_loc + tral(:, ((i-1)*i + 1):end);     % form the right-hand side i copies
    paraupp.r = src_loc + flipud(tral(:, ((i-1)*i + 1):end));
    paraded.r = src_loc + flipud(tral(:,1:i));
    
    paralef.n = [1;0];
    pararig.n = [1;0];
    paraupp.n = [1;0];
    paraded.n = [1;0];
    
    e_2 = [0;1];
    e_1 = [1;0];

    Ded = chnkr1; Upp = chnkr1.move([], e_2); 

    Lef = chnkr2; Rig = chnkr2.move([], e_1);


    valslefsrc = 0;                                             % first discrepancy condition
    for k = 1:length(paralef(1,:))
    valslefsrc = valslefsrc + sum(kern1(paralef(:, k), Rig),2);
    end

    valsrigsrc = 0;
    for k = 1:length(pararig(1,:))
    valsrigsrc = valsrigsrc + sum(kern1(pararig(:, k), Lef),2);
    end

    vals2lefsrc = 0;                                            % second discrepancy condition
    for k = 1:length(paralef(1,:))
    vals2lefsrc = vals2lefsrc + sum(kern2( paralef(:, k), Rig ),2);
    end
    
    vals2rigsrc = 0;
    for k = 1:length(pararig(1,:))
    vals2rigsrc = vals2rigsrc + sum(kern2( pararig(:,k), Lef),2);
    end

    valsdedsrc = 0;                                             % third discrepancy condition
    for k = 1:length(paraded(1,:))
    valsdedsrc= valsdedsrc + sum(kern1( paraded(:,k), Upp),2);
    end

    valsuppsrc = 0;
    for k = 1:length(paraupp(1,:))
    valsuppsrc  = valsuppsrc  + sum(kern1( paraupp(:,k), Ded),2);
    end

    vals2dedsrc = 0;                                             % forth discrepancy condition
    for k = 1:length(paraded(1,:))
    vals2dedsrc= vals2dedsrc + sum(kern2(paraded(:, k), Upp ),2);
    end

    vals2uppsrc = 0;
    for k = 1:length(paraupp(1,:))
    vals2uppsrc  = vals2uppsrc  + sum(kern2(paraupp(:,k) ,Ded),2);
    end

    C1 = valslefsrc - valsrigsrc; 
    C2 = vals2lefsrc - vals2rigsrc;             % note the sign flip due to the different normal convention
    C3 = valsdedsrc - valsuppsrc;
    C4 = vals2dedsrc - vals2uppsrc;

    Cblockrhs = [C1; C2; C3; C4];
end