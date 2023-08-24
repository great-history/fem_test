%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 得到标准一维[-1,1]的节点上的形函数的偏导
% get the standard shape functions' partial x/y matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function std_shape_list = get_std_1d_shape_mat(r)
    std_shape_list = zeros(2, 1);
    std_shape_list(1, 1) = 0.5 * (1 - r);
    std_shape_list(2, 1) = 0.5 * (1 + r);
end