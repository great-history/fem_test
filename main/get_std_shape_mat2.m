%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 得到标准三角形的节点上的形函数的值
% get the standard shape functions' value matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function std_shape_mat = get_std_shape_mat2(r, s, num_node_per_unit)
    switch num_node_per_unit
        case 3
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 对于平面3节点三角形情形
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            std_shape_mat = zeros(1, 3);
            
            std_shape_mat(1, 1) = 1 - r - s;
            std_shape_mat(1, 2) = r;
            std_shape_mat(1, 3) = s;
        case 6
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 对于平面6节点三角形情形 (待补充)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            std_shape_mat = zeros(1, 6);
            
            std_shape_mat(1, 1) = (1 - 2 * r - 2 * s) * (1 - r - s);
            std_shape_mat(1, 2) = (2 * r - 1) * r;
            std_shape_mat(1, 3) = (2 * s - 1) * s;
            std_shape_mat(1, 4) = 4 * (1 - r - s) * r;
            std_shape_mat(1, 5) = 4 * r * s;
            std_shape_mat(1, 6) = 4 * s * (1 - r - s);
        otherwise
            error('too many nodes per triangular unit !!!')
    end
end