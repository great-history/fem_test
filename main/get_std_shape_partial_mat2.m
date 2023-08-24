%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 得到标准三角形的节点上的形函数的偏导
% get the standard shape functions' partial x/y matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function std_shape_partial_mat = get_std_shape_partial_mat2(r, s, num_node_per_unit)
    switch num_node_per_unit
        case 3
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 对于平面3节点三角形情形
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            std_shape_partial_mat = zeros(2, 3);
            std_shape_partial_mat(1, 1) = - 1;
            std_shape_partial_mat(2, 1) = - 1;

            std_shape_partial_mat(1, 2) = 1;
            std_shape_partial_mat(2, 2) = 0;

            std_shape_partial_mat(1, 3) = 0;
            std_shape_partial_mat(2, 3) = 1;
        case 6
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 对于平面6节点三角形情形 (待补充)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            std_shape_partial_mat = zeros(2, 6);
            
            std_shape_partial_mat(1, 1) = - 3 + 4 * r + 4 * s;
            std_shape_partial_mat(2, 1) = - 3 + 4 * r + 4 * s;
            
            std_shape_partial_mat(1, 2) = 4 * r - 1;
            std_shape_partial_mat(2, 2) = 0;
            
            std_shape_partial_mat(1, 3) = 0;
            std_shape_partial_mat(2, 3) = 4 * s - 1;
            
            std_shape_partial_mat(1, 4) = 4 * (1 - 2 * r - s);
            std_shape_partial_mat(2, 4) = - 4 * r;
            
            std_shape_partial_mat(1, 5) = 4 * s;
            std_shape_partial_mat(2, 5) = 4 * r;
            
            std_shape_partial_mat(1, 6) = - 4 * s;
            std_shape_partial_mat(2, 6) = 4 * (1 - r - 2 * s);
        otherwise
            error('too many nodes per triangular unit !!!')
    end
end