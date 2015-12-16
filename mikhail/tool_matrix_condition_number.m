%> @file tool_matrix_condition_number.m
%> @brief This function calculates the Matrix Condition Number.
%> @author Mikhail V. Konnik
%> @date   25 May 2011
%>
%> @section condnumebr Matrix Condition Number
%> By the definition of the singular values of a matrix, An [n x n] matrix P has n singular values \f$\sigma\f$, given as \f$\sigma^2 = \lambda(P^T P) = \lambda(P P^T)\f$. The matrix \f$(P^TP) is symmetric, and the eigenvalues of a symmetric matrix are always real and non-negative, so the singular values of a matrix are always real and non-negative. The matrix P is non-singular (invertible) if and only if all of its singular values are positive.
%>
%> The condition number of a matrix is defined as:@n
%> \f$\kappa(P) = \frac{\sigma_{max}(P)}{\sigma_{min}(P)}\f$
%>
%> As \f$\kappa(P) \rightarrow \infty\f$, the matrix P is said to be poorly conditioned or ill conditioned, and P approaches a singular matrix.
%======================================================================
%> @param A 	= the matrix whose condition number must be calculated.
%> @retval condition_number = the condition number for the input matrix A.
% ======================================================================
function [condition_number,S1] = tool_matrix_condition_number(A);

if ( (min(size(A)))<= 2)
	disp('The dimension of matrix must be greater than 2.');
	condition_number = NaN;
	return;
end %% if (min(size(A)))<= 2)


if (issparse(A) == 1)
	[U,S,V] = svds(A); %% performing the Singular Values Decomposition (SVD).
else
	[U,S,V] = svd(A); %% performing the Singular Values Decomposition (SVD).
end %% if (issparse(A)) ==1)

S1 = diag(S); %%% extracting only the diagonal elements from the S matrix of singular values.

condition_number = max(S1)/min(S1);