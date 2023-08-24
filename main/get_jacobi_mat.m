%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 得到对应(r,s)的Jacobi矩阵
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function jacobi_mat = get_jacobi_mat(r, s, unit_node_pos_list)
    std_shape_partial_mat = get_std_shape_partial_mat(r, s);
    jacobi_mat = std_shape_partial_mat * unit_node_pos_list;
end