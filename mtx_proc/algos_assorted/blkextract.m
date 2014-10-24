%**************************************************************************
% Extracts unique blocks (no overlap) of size mxn from the input matrix.
%
% @param   A   The matrix to extract the data from.
% @param   m   The block row dimension.
% @param   n   The block column dimension.
%
% @return  A cell matrix of blocks.
%
% @note    If any block falls outside of the matrix, then the block is
%          padded with zeros in that area that falls outside.
%
% @throw   Error if there are more than two dimensions in the matrix.
%
%**************************************************************************
% Code is written by author based upon noted references (if given), using 
% The MathWorks MATLAB function signature (if applicable) for convenience 
% only.
%
% If you have any questions, comments, or find bugs, please feel free to 
% email me at geoff.hayes74@gmail.com.
%
% Geoff Hayes 2014
%**************************************************************************
function [blks] = blkextract(A,m,n)

    if ~ismatrix(A)
        error('blkextract: more than two dimensions is unsupported');
    end
    
    mb = ceil(abs(m(1)));
    nb = ceil(abs(n(1)));

    [ma,na] = size(A);  % get the number of rows and cols in x

    rBlks = floor((ma-1)/mb) + 1;  % determine the number of row blocks
    cBlks = floor((na-1)/nb) + 1;  % determine the number of col blocks

    blks = cell(rBlks,cBlks);  % size your output cell matrix

    % extract all blocks from the matrix x
    for i=1:rBlks

        % determine the start row (xra) and the end row (xrb)
        % for row block i
        xra = (i-1)*mb+1; xrb = min(xra+mb-1,ma);

        for j=1:cBlks

            % determine the start col (xca) and the end col (xcb)
            % for column block j
            xca = (j-1)*nb+1; xcb = min(xca+nb-1,na);

            % reset the empty block to all zeros (zeros will be the
            % padding in the event the block we extract from x is
            % not complete)
            blk = zeros(mb,nb);

            % extract the data (which may now be padded)
            blk(1:(xrb-xra+1),1:(xcb-xca+1)) = A(xra:xrb,xca:xcb);

            % save to cell matrix
            blks{i,j} = blk;
        end
    end
end
