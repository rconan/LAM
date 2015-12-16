function tool_inside_path_adder(this_project_name, other_project_name)

current_full_path = pwd; %%% get the full path to the currently running function
general_path = strrep(current_full_path, this_project_name, '');

addpath(strcat(general_path, other_project_name) );