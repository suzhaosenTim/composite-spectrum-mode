function [objcchnkr,Cblock] = formCblockmat(chnkobj, ctr, kern1, kern2, chnkr1, chnkr2)
% A method to form the C Block matrix;
% Input argument: chnkobj- a chunker or chunkgraph describing the inclusion
% ctr: center of the object
% kern1: first kernel- typically G(x,y) 
% kern2: second kernel- typically \partial_n_x G(x,y)
% 
% Note: self-evaluation is handled separately. 


if (class(chnkobj) == "chunker" )
    chnkr = chnkobj;
    [x, y] = meshgrid(-1:1); tral = [x(:).'; y(:).']; 

                    
     
    paralef = ctr + tral(:,1:3);       % form the left-hand side 3 copies
    pararig = ctr + tral(:, 7:end);     % form the right-hand side 3 copies
    paraupp = ctr + tral(:, [3, 6, 9]);
    paraded = ctr + tral(:, [1, 4, 7]);

    objclef = sumSrccopies(chnkr, paralef);
    
    
    objcrig = sumSrccopies(chnkr, pararig);     % build left and right chnkr
    
    objcupp = sumSrccopies(chnkr, paraupp);
    
    objcded = sumSrccopies(chnkr, paraded);    % build up and down chnkr
    
    % used to check what Srccopies built us.
    objcchnkr = merge([objclef, objcrig, objcupp, objcded]);
    
    e_2 = [0;1];
    e_1 = [1;0];

    Ded = chnkr1; Upp = chnkr1.move([], e_2); 

    Lef = chnkr2; Rig = chnkr2.move([], e_1);
    


    valslefsrc = 0;                                             % first discrepancy condition
    for k = 1:numel(objclef)
    weights1 = objclef(k).wts(:).';
    valslefsrc = valslefsrc + kern1(objclef(k), Rig).*weights1;
    end

    valsrigsrc = 0;
    for k = 1:numel(objcrig)
    weights2 = objcrig(k).wts(:).';
    valsrigsrc = valsrigsrc + kern1(objcrig(k), Lef).*weights2;
    end

    vals2lefsrc = 0;                                            % second discrepancy condition
    for k = 1:numel(objclef)
    weights3 = objclef(k).wts(:).';
    vals2lefsrc = vals2lefsrc + kern2(objclef(k), Rig).*weights3;
    end
    
    vals2rigsrc = 0;
    for k = 1:numel(objcrig)
    weights4 = objcrig(k).wts(:).';
    vals2rigsrc = vals2rigsrc + kern2(objcrig(k), Lef).*weights4;
    end

    valsdedsrc = 0;                                             % third discrepancy condition
    for k = 1:numel(objcded)
    weights5 = objcded(k).wts(:).';
    valsdedsrc= valsdedsrc + kern1(objcded(k), Upp).*weights5;
    end

    valsuppsrc = 0;
    for k = 1:numel(objcupp)
    weights6 = objcupp(k).wts(:).';
    valsuppsrc  = valsuppsrc  + kern1(objcupp(k), Ded).*weights6;
    end

    vals2dedsrc = 0;                                             % forth discrepancy condition
    for k = 1:numel(objcded)
    weights7 = objcded(k).wts(:).';
    vals2dedsrc= vals2dedsrc + kern2(objcded(k), Upp).*weights7;
    end

    vals2uppsrc = 0;
    for k = 1:numel(objcupp)
    weights8 = objcupp(k).wts(:).';
    vals2uppsrc  = vals2uppsrc  + kern2(objcupp(k), Ded).*weights8;
    end

    C1 = valslefsrc - valsrigsrc; 
    C2 = vals2lefsrc - vals2rigsrc;             % note the sign flip due to the different normal convention
    C3 = valsdedsrc - valsuppsrc;
    C4 = vals2dedsrc - vals2uppsrc;

    Cblock = [C1; C2; C3; C4];

end
    
    
end