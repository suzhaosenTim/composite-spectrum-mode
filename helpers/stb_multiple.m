function mat =  stb_multiple(chnkr_int)
% helper method to build stabilizer for multiply-connected geometry. 
% Helpful for the Laplace Neumann problem. 
% input: an array of chunker objects (default chnkr, needs to be tested for
% cgrphs)
% output: blk-diag matrix with \int_{\gamma_{i}}\mu dS on the diagonal. 

stb= {};
for i = 1:length(chnkr_int)
    stb{i} = onesmat(chnkr_int(i));
end

mat = blkdiag(stb{:});

end
