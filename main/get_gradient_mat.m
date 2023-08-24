%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 得到待求函数的偏导（数值）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gradient_mat = get_gradient_mat(r, s, unit_node_value_list)
    std_shape_partial_mat = get_std_shape_partial_mat(r, s);
    gradient_mat = std_shape_partial_mat * unit_node_value_list;
end