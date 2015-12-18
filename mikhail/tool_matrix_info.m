%> @file tool_matrix_info.m
%> @brief The tool for the analysis of the properties of the matrix.
%>
%> @author Mikhail Konnik
%> @date   7 August 2012
%>
%> @section matpropers Matrix Properties
%> It turns out that the way how the properties of the matrix are stated in textbooks are very different from how the properties (det, eig...) are @b actually calculated.
%> This function uses both approaches and returns the structure of the values of many matrix properties, including condition number, determinant, eigs, svd and all the rest.
%======================================================================
%> @param structure_to_print        = the structure to be printed into the file.
%> @param name_to_print             = output file name for the structure.
% ======================================================================
function matrix = tool_matrix_info(varargin);

switch nargin
    case 1
        Q = varargin{1};
    case 2
        Q = varargin{1};
        name_to_print = varargin{2};
    otherwise
        error('Error in tool_matrix_info: first arg must be a matrix, and the second - the filename for output (optional).')
end 

warning off

[matrix.rows,  matrix.cols] = size(Q);


%%%%%%%% Real/Complex matrix?
if (nnz(imag(Q)) ~= 0)
    matrix.type = 'Complex';
else
    matrix.type = 'Real';
end


%%%%%%%% Class Sparse/Dense matrix?
if (issparse(Q) ~= 0)
    matrix.class = 'Sparse';
else
    matrix.class = 'Dense';
end


%%%% Call tool for the symmetry %%%
if (matrix.rows == matrix.cols)
    matrix.symmetry = tool_issymmetric(Q);

    if matrix.symmetry == 0
        matrix.symmetry = 'Non-Symmetric'
    else
        if ( matrix.type == 'Real')
            matrix.symmetry = 'Symmetric';
        else
            matrix.symmetry = 'Hermitian';
        end
    end
else
    matrix.symmetry = 'Non-Symmetric (Rectangular)';
end




%%%% Basic properties %%%
if strcmp(matrix.symmetry , 'Non-Symmetric (Rectangular)')==1;
    matrix.determinant = 'N/A (Matrix must be square)';
else
    matrix.determinant = det(Q);
%     
% Algorithm
% 
% 
% The determinant is computed from the triangular factors obtained by Gaussian elimination
% [L,U] = lu(A)
% s =  det(L)        % This is always +1 or -1 
% det(A) = s*prod(diag(U))
end

% matrix.trace = trace(Q);


%%%%% Sparsity %%%%%
input_matrix = sparse(Q); %%% force the matrix to be sparse

matrix.total_elements = size(input_matrix,1)*size(input_matrix,2);
matrix.non_zero_elements = nnz(input_matrix(:));
matrix.sparsity_percent = round( ((matrix.total_elements - matrix.non_zero_elements)/matrix.total_elements )*100) ;


% 
% 
% Some operations on a vector
% x = [11  0  33  0  55]';
% 
% find(x)
% ans =
%      1
%      3
%      5
% 
% find(x == 0)
% ans =
%      2
%      4
% 
% find(0 < x & x < 10*pi)
% ans =
%      1


if strcmp(matrix.symmetry , 'Non-Symmetric (Rectangular)')==1;
    matrix.eig = 'N/A (Matrix must be square)';
else
    
%%%%%% Eigevalues %%%%
% if issparse(Q) == 1
    [U,V] = eig(full(Q));
% else
% %     [U,V] = eig(Q);
% end

matrix.eig_min = min(diag(V));
matrix.eig_max = max(diag(V));

matrix.spectral_radius = abs( max(diag(V)));

matrix.eig_tolernace_eps = 1;
matrix.eig_negatives = nnz(find(diag(V)<= -matrix.eig_tolernace_eps*eps));
matrix.eig_epssize = nnz(find(-matrix.eig_tolernace_eps*eps <= diag(V) <= matrix.eig_tolernace_eps*eps));

matrix.eig_positives = nnz(find(diag(V)> matrix.eig_tolernace_eps*eps));


% matrix.eig_condition_number_condeig = condeig(full(Q));


%%%%%%%%%%%%% Positive Definiteness? %%%%%
if strcmp(matrix.symmetry,'Non-Symmetric') ==1;
    matrix.definiteness = 'Indefinite';
else
    if (matrix.eig_negatives>0)
        matrix.definiteness = 'Indefinite';
    elseif (matrix.eig_negatives==0 && matrix.eig_epssize>0)
        matrix.definiteness = 'Positive Semi-definite';
    elseif (matrix.eig_negatives==0 && matrix.eig_epssize==0 && matrix.eig_positives>0)
        matrix.definiteness = 'Positive Definite';
    end
end


end %%% Is the matrix is square



%%%%%%%%%%%%% Condition number and SVD %%%%%
[condition_number,S1] = tool_matrix_condition_number(Q);

matrix.singular_val_min = min(S1);
matrix.singular_val_max = max(S1);

 matrix.condition_number_normwise2 = matrix.singular_val_max/matrix.singular_val_min;

%%%%%%%%%%%%% Condition number estimatoin %%%%%
if (matrix.rows == matrix.cols)
    matrix.condition_number_matlab_norm2 = cond(Q,2);
    matrix.condition_number_matlab_condest = condest(Q);
    matrix.condition_number_matlab_rcon = rcond(full(Q));
else
    matrix.condition_number_matlab = 'N/A, Matrix must be square';
end



%%%%%%%%%% Rank of the matrix %%%%%%%%%%
% There are a number of ways to compute the rank of a matrix. 
% MATLAB uses the method based on the singular value decomposition, or SVD. 
% The SVD algorithm is the most time consuming, but also the most reliable.
tol = max(size(Q))*eps(max(S1));
r = sum(S1 > tol);
matrix.rank = r;






%%%%%%%%%%%%%%%%%%%%%%% If a user wants to print all those properties out - let it be %%%%%%%%%%%%%%
if nargin == 2
    tool_print_structure_contents (matrix,name_to_print);
end
%%%%%%%%%%%%%%%%%%%%%%% If a user wants to print all those properties out - let it be %%%%%%%%%%%%%%