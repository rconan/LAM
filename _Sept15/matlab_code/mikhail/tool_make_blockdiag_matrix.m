function ZBlockDiag = tool_make_blockdiag_matrix (matrix_name_string, nGs, matrix)

eval_command_tmp = matrix_name_string;

for ii = 2:nGs
    eval_command_tmp = strcat(eval_command_tmp, ',' , matrix_name_string);
end

eval_temp_string = strcat (matrix_name_string, '= matrix;');
eval(eval_temp_string); %% evaluate the command as a string


eval_command = strcat('ZBlockDiag  = blkdiag(', eval_command_tmp, ');');
eval(eval_command); %% evaluate the command as a string

% % % CnAst = blkdiag( Cn, Cn, Cn ); %%% block-diagonal Covariance matrix, the size of N Ngs