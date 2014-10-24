% Function to write binary data files of random matrices in order to test
% the TMatrix<class T> C++ matrix library.  Note that the transpose of each
% matrix is written to file.  This is because the MATLAB fwrite writes out
% the matrix data in columnn order, whereas the C++ file reader expects the
% data by rows.
%
% Input:
%   saveFilesToDir   The directory to save the binary data files to.
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
function tmatrix_driver(saveFilesToDir)

    % reset the random number generator
    rng('default');

    % close all open files
    fclose all;
    
    % maximum dimension for a matrix row or column
    maxMatrixDim         = 10;
    
    % maximum number of matrices to generate per (row,column) pair
    maxNumMatricesPerDim = 100;

    % define data precision: single or double
    doSinglePrecision=0;
    doDoublePrecision=1;
    
    if doSinglePrecision
        dataTypeDir = 'float';
    else
        dataTypeDir = 'double';
    end
    
    saveFilesToDir = [saveFilesToDir '/' dataTypeDir '/'];
    
    if ~isdir(saveFilesToDir)
        fprintf('The directory %s is invalid!!\n', saveFilesToDir);
        return;
    end   
    
    % test out matrix assignment and the copy constructor
    fod = fopen([saveFilesToDir 'mtxAssignmentTest.bin'], 'wb');
    
    if fod
        
        % generate randomly sized matrices (varying the rows and columns)
        % up to the max matrix dimension
        for i=1:maxMatrixDim
            for j=1:maxMatrixDim
                for k=1:maxNumMatricesPerDim
                   
                    % generate a matrix given the dimensions
                    A = generateMatrix(i,j,doSinglePrecision);
                    
                    % write the matrix dimensions and the matrix to file
                    fwrite(fod, i, 'uint16');
                    fwrite(fod, j, 'uint16');
                    if doSinglePrecision
                        fwrite(fod,A,'float');
                    else
                        fwrite(fod,A,'double');
                    end
                
                end
            end
        end
        
        fclose(fod);
    end
    
    % test out matrix scalar assignment
    fod = fopen([saveFilesToDir 'mtxScalarAssignmentTest.bin'], 'wb');
    
    if fod
        
        % generate randomly sized matrices (varying the rows and columns)
        % up to the max matrix dimension
        for i=1:maxMatrixDim
            for j=1:maxMatrixDim
                for k=1:maxNumMatricesPerDim
                   
                    % generate a scalar assignment value
                    scalar = generateMatrix(1,1,doSinglePrecision);
                    
                    % write the matrix dimensions and the 1x1 "matrix" to 
                    % file
                    fwrite(fod, i, 'uint16');
                    fwrite(fod, j, 'uint16');
                    if doSinglePrecision
                        fwrite(fod,scalar,'float');
                    else
                        fwrite(fod,scalar,'double');
                    end
                
                end
            end
        end
        
        fclose(fod);
    end

    % test out matrix scalar addition
    fod = fopen([saveFilesToDir 'mtxScalarAdditionTest.bin'], 'wb');
    
    if fod
        
        % generate randomly sized matrices (varying the rows and columns)
        % up to the max matrix dimension
        for i=1:maxMatrixDim
            for j=1:maxMatrixDim
                for k=1:maxNumMatricesPerDim
                   
                    % generate a scalar assignment value
                    scalar = generateMatrix(1,1,doSinglePrecision);
                    
                    % generate a matrix
                    A = generateMatrix(i,j,doSinglePrecision);
                    
                    % add the scalar to the matrix
                    B = A + scalar;
                    
                    % write the data to file
                    fwrite(fod, i, 'uint16');
                    fwrite(fod, j, 'uint16');
                    if doSinglePrecision
                        fwrite(fod,scalar,'float');
                        fwrite(fod,A','float');
                        fwrite(fod,B','float');
                    else
                        fwrite(fod,scalar,'double');
                        fwrite(fod,A','double');
                        fwrite(fod,B','double');
                    end
                end
            end
        end
        
        fclose(fod);
    end
    
    % test out matrix scalar subtraction
    fod = fopen([saveFilesToDir 'mtxScalarSubtractionTest.bin'], 'wb');
    
    if fod
        
        % generate randomly sized matrices (varying the rows and columns)
        % up to the max matrix dimension
        for i=1:maxMatrixDim
            for j=1:maxMatrixDim
                for k=1:maxNumMatricesPerDim
                   
                    % generate a scalar assignment value
                    scalar = generateMatrix(1,1,doSinglePrecision);
                    
                    % generate a matrix
                    A = generateMatrix(i,j,doSinglePrecision);
                    
                    % add the scalar to the matrix
                    B = A - scalar;
                    
                    % write the data to file
                    fwrite(fod, i, 'uint16');
                    fwrite(fod, j, 'uint16');
                    if doSinglePrecision
                        fwrite(fod,scalar,'float');
                        fwrite(fod,A','float');
                        fwrite(fod,B','float');
                    else
                        fwrite(fod,scalar,'double');
                        fwrite(fod,A','double');
                        fwrite(fod,B','double');
                    end
                end
            end
        end
        
        fclose(fod);
    end    
    
    % test out matrix scalar multiplication
    fod = fopen([saveFilesToDir 'mtxScalarMultiplicationTest.bin'], 'wb');
    
    if fod
        
        % generate randomly sized matrices (varying the rows and columns)
        % up to the max matrix dimension
        for i=1:maxMatrixDim
            for j=1:maxMatrixDim
                for k=1:maxNumMatricesPerDim
                   
                    % generate a scalar assignment value
                    scalar = generateMatrix(1,1,doSinglePrecision);
                    
                    % generate a matrix
                    A = generateMatrix(i,j,doSinglePrecision);
                    
                    % add the scalar to the matrix
                    B = A * scalar;
                    
                    % write the data to file
                    fwrite(fod, i, 'uint16');
                    fwrite(fod, j, 'uint16');
                    if doSinglePrecision
                        fwrite(fod,scalar,'float');
                        fwrite(fod,A','float');
                        fwrite(fod,B','float');
                    else
                        fwrite(fod,scalar,'double');
                        fwrite(fod,A','double');
                        fwrite(fod,B','double');
                    end
                end
            end
        end
        
        fclose(fod);
    end
    
    % test out matrix scalar division
    fod = fopen([saveFilesToDir 'mtxScalarDivisionTest.bin'], 'wb');
    
    if fod
        
        % generate randomly sized matrices (varying the rows and columns)
        % up to the max matrix dimension
        for i=1:maxMatrixDim
            for j=1:maxMatrixDim
                for k=1:maxNumMatricesPerDim
                   
                    % generate a scalar assignment value
                    scalar = generateMatrix(1,1,doSinglePrecision);
                    
                    % generate a matrix
                    A = generateMatrix(i,j,doSinglePrecision);
                    
                    % add the scalar to the matrix
                    B = A / scalar;
                    
                    % write the data to file
                    fwrite(fod, i, 'uint16');
                    fwrite(fod, j, 'uint16');
                    if doSinglePrecision
                        fwrite(fod,scalar,'float');
                        fwrite(fod,A','float');
                        fwrite(fod,B','float');
                    else
                        fwrite(fod,scalar,'double');
                        fwrite(fod,A','double');
                        fwrite(fod,B','double');
                    end
                end
            end
        end
        
        fclose(fod);
    end
    
    % test out matrix addition
    fod = fopen([saveFilesToDir 'mtxAdditionTest.bin'], 'wb');
    
    if fod
        
        % generate randomly sized matrices (varying the rows and columns)
        % up to the max matrix dimension
        for i=1:maxMatrixDim
            for j=1:maxMatrixDim
                for k=1:maxNumMatricesPerDim
                   
                    % generate two matrices
                    A = generateMatrix(i,j,doSinglePrecision);
                    B = generateMatrix(i,j,doSinglePrecision);
                    
                    % add the two matrices
                    C = A + B;
                    
                    % write the data to file
                    fwrite(fod, i, 'uint16');
                    fwrite(fod, j, 'uint16');
                    if doSinglePrecision
                        fwrite(fod,A','float');
                        fwrite(fod,B','float');
                        fwrite(fod,C','float');
                    else
                        fwrite(fod,A','double');
                        fwrite(fod,B','double');
                        fwrite(fod,C','double');
                    end
                end
            end
        end
        
        fclose(fod);
    end    

    % test out matrix subtraction
    fod = fopen([saveFilesToDir 'mtxSubtractionTest.bin'], 'wb');
    
    if fod
        
        % generate randomly sized matrices (varying the rows and columns)
        % up to the max matrix dimension
        for i=1:maxMatrixDim
            for j=1:maxMatrixDim
                for k=1:maxNumMatricesPerDim
                   
                    % generate two matrices
                    A = generateMatrix(i,j,doSinglePrecision);
                    B = generateMatrix(i,j,doSinglePrecision);
                    
                    % subtract the two matrices
                    C = A - B;
                    
                    % write the data to file
                    fwrite(fod, i, 'uint16');
                    fwrite(fod, j, 'uint16');
                    if doSinglePrecision
                        fwrite(fod,A','float');
                        fwrite(fod,B','float');
                        fwrite(fod,C','float');
                    else
                        fwrite(fod,A','double');
                        fwrite(fod,B','double');
                        fwrite(fod,C','double');
                    end
                end
            end
        end
        
        fclose(fod);
    end      

    % test out matrix multiplication
    fod = fopen([saveFilesToDir 'mtxMultiplicationTest.bin'], 'wb');
    
    if fod
        ctr = 0;
        % generate randomly sized matrices (varying the rows and columns)
        % up to the max matrix dimension
        for i=1:maxMatrixDim
            for j=1:maxMatrixDim
                for k=1:maxNumMatricesPerDim

                    % generate a random column dimension for the second
                    % matrix
                    m = 1 + floor(rand(1)*maxMatrixDim);
                   
                    % generate two matrices
                    A = generateMatrix(i,j,doSinglePrecision);
                    B = generateMatrix(j,m,doSinglePrecision);
                    
                    % multiply the two matrices
                    C = A * B;
                    
                    % write the data to file
                    fwrite(fod, i, 'uint16');
                    fwrite(fod, j, 'uint16');
                    fwrite(fod, j, 'uint16');
                    fwrite(fod, m, 'uint16');
                                        
                    if doSinglePrecision
                        fwrite(fod,A','float');
                        fwrite(fod,B','float');
                        fwrite(fod,C','float');
                    else
                        fwrite(fod,A','double');
                        fwrite(fod,B','double');
                        fwrite(fod,C','double');
                    end
                end
            end
        end
        
        fclose(fod);
    end  
    
    % test out matrix transpose
    fod = fopen([saveFilesToDir 'mtxTransposeTest.bin'], 'wb');
    
    if fod
        ctr = 0;
        % generate randomly sized matrices (varying the rows and columns)
        % up to the max matrix dimension
        for i=1:maxMatrixDim
            for j=1:maxMatrixDim
                for k=1:maxNumMatricesPerDim
                  
                    % generate the matrix
                    A = generateMatrix(i,j,doSinglePrecision);
                    
                    % write the data to file
                    fwrite(fod, i, 'uint16');
                    fwrite(fod, j, 'uint16');
                                
                    % write out the matrix and its transpose
                    if doSinglePrecision
                        fwrite(fod,A','float');
                        fwrite(fod,A,'float');
                    else
                        fwrite(fod,A','double');
                        fwrite(fod,A,'double');
                    end
                end
            end
        end
        
        fclose(fod);
    end  

    % test out matrix determinant
    fod = fopen([saveFilesToDir 'mtxDeterminantTest.bin'], 'wb');
    
    if fod
        ctr = 0;
        % generate randomly sized square matrices (varying the rows and
        % columns) up to the max matrix dimension
        for i=1:maxMatrixDim
            for j=1:maxMatrixDim
                for k=1:maxNumMatricesPerDim
                  
                    % generate the square matrix
                    A = generateMatrix(i,i,doSinglePrecision);
                    
                    % write the data to file
                    fwrite(fod, i, 'uint16');
                    fwrite(fod, i, 'uint16');
                                
                    % write out the matrix and its determinant
                    if doSinglePrecision
                        fwrite(fod,A','float');
                        fwrite(fod,det(A),'float');
                    else
                        fwrite(fod,A','double');
                        fwrite(fod,det(A),'double');
                    end
                end
            end
        end
        
        fclose(fod);
    end  
    
    % test out matrix inverse
    fod = fopen([saveFilesToDir 'mtxInverseTest.bin'], 'wb');
    
    if fod
        ctr = 0;
        % generate randomly sized square matrices (varying the rows and
        % columns) up to the max matrix dimension
        for i=1:maxMatrixDim
            for j=1:maxMatrixDim
                for k=1:maxNumMatricesPerDim
                  
                    % generate the square matrix
                    A = generateMatrix(i,i,doSinglePrecision);
                    
                    while det(A)==0
                        A = generateMatrix(i,i,doSinglePrecision);
                    end                        
                    
                    % write the data to file
                    fwrite(fod, i, 'uint16');
                    fwrite(fod, i, 'uint16');
                                
                    % write out the matrix and its determinant
                    if doSinglePrecision
                        fwrite(fod,A','float');
                        fwrite(fod,inv(A)','float');
                    else
                        fwrite(fod,A','double');
                        fwrite(fod,inv(A)','double');
                    end
                end
            end
        end
        
        fclose(fod);
    end  
    
    % test out matrix LUP
    fod = fopen([saveFilesToDir 'mtxLUPTest.bin'], 'wb');
    
    if fod
        ctr = 0;
        % generate randomly sized square matrices (varying the rows and
        % columns) up to the max matrix dimension
        for i=1:maxMatrixDim
            for j=1:maxMatrixDim
                for k=1:maxNumMatricesPerDim
                  
                    % generate the square matrix
                    A = generateMatrix(i,i,doSinglePrecision);
                    
                    % calculate the LUP
                    [L,U,P] = lu(A);

                    % write the data to file
                    fwrite(fod, i, 'uint16');
                    fwrite(fod, i, 'uint16');
                                
                    % write out the matrix and its determinant
                    if doSinglePrecision
                        fwrite(fod,A','float');
                        fwrite(fod,L','float');
                        fwrite(fod,U','float');
                        fwrite(fod,P','float');
                    else
                        fwrite(fod,A','double');
                        fwrite(fod,L','double');
                        fwrite(fod,U','double');
                        fwrite(fod,P','double');
                    end
                end
            end
        end
        
        fclose(fod);
    end     
end

% Function to generate random entries within a matrix given the required
% number of rows and columsn
function [matrix] = generateMatrix(numRows, numCols, doSinglePrecision)

    if doSinglePrecision
        a = single(-100);
        b = single(100);
        matrix = a + (b-a).*single(rand(numRows,numCols));
    else
        a = -100;
        b = 100;
        matrix = a + (b-a).*rand(numRows,numCols);
    end
    
end


