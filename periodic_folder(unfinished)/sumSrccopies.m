function objc = sumSrccopies(chnkobj, tral)
% A method to create chnkr/chunkgraphs n nearby source copies
% and return an array contains these objects


if (class(chnkobj) == "chunker" )
    
    chnkr = chnkobj;
    objc = [];
    for k = 1:length(tral)
        objc = [objc, chnkr.move([], tral(:,k))];
    end
end

if (class(chnkobj) == "chunkgraph" )        % This needs to be tested later. 
    
    cgrph = chnkobj;
    objc = [];
    for k = 1:length(tral)
        objc = [objc, cgrph.move([], tral(:,k))];
    end
end



end