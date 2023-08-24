%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 得到标准四边形的节点上的形函数的偏导
% get the standard shape functions' partial x/y matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function std_shape_partial_mat = get_std_shape_partial_mat(r, s)
    std_shape_partial_mat = zeros(2, 4);
    std_shape_partial_mat(1, 1) = - 0.25 * (1 - s);
    std_shape_partial_mat(2, 1) = - 0.25 * (1 - r);
    
    std_shape_partial_mat(1, 2) =   0.25 * (1 - s);
    std_shape_partial_mat(2, 2) = - 0.25 * (1 + r);
    
    std_shape_partial_mat(1, 3) =   0.25 * (1 + s);
    std_shape_partial_mat(2, 3) =   0.25 * (1 + r);
    
    std_shape_partial_mat(1, 4) = - 0.25 * (1 + s);
    std_shape_partial_mat(2, 4) =   0.25 * (1 - r);
end