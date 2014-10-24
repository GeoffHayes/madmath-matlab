%**************************************************************************
% Extracts overlapping blocks of size mxm from the input matrix.
%
% @param   A   The matrix to extract the data from.
% @param   m   The square block dimension.
% @param   s   The block centre shift parameter (if s==m, then there is no
%              overlap).
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
function [blks] = ovlpblkextract(A,m,s)

    if ndims(A)>3
        error('ovlpblkextract: more than three dimensions is unsupported');
    end
    
    ctype = class(A);
    
    mb = ceil(abs(m(1)));
    
    doOverlap = false;
    
    if ~exist('s','var')
        % no overlap
        s = mb;
    else
        % s can never be greater than the block width else we would be
        % missing elements from A
        s = min(mb,ceil(abs(s(1))));
        
        if s~=mb
            % there is overlap, ensure that mb is odd
            if mod(mb,2)==0
                mb = mb + 1;
            end
        end
        
        doOverlap = true;
    end
    
    [ma,na,pa] = size(A);  % get the number of rows and cols in x
    
    if ~doOverlap
        rBlks = floor((ma-1)/s) + 1;  % determine the number of row blocks
        cBlks = floor((na-1)/s) + 1;  % determine the number of col blocks
    else
        rBlks = floor((ma-1-(ceil(mb/2)))/s) + 1;  % determine the number of row blocks
        cBlks = floor((na-1-(ceil(mb/2)))/s) + 1;  % determine the number of col blocks
    end

    blks = cell(rBlks,cBlks);  % size your output cell matrix

    % extract all blocks from the matrix x
    for i=1:rBlks

        % determine the start row (xra) and the end row (xrb)
        % for row block i
        xra = (i-1)*s+1; xrb = min(xra+mb-1,ma);

        for j=1:cBlks

            % determine the start col (xca) and the end col (xcb)
            % for column block j
            xca = (j-1)*s+1; xcb = min(xca+mb-1,na);

            % reset the empty block to all zeros (zeros will be the
            % padding in the event the block we extract from x is
            % not complete)
            blk = cast(zeros(mb,mb,pa),ctype);

            % extract the data (which may now be padded)
            blk(1:(xrb-xra+1),1:(xcb-xca+1),:) = A(xra:xrb,xca:xcb,:);

            % save to cell matrix
            blks{i,j} = blk;
        end
    end
end
