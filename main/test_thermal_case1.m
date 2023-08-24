%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 读取网格(这里采用三角形划分)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
node_list = h5read('thermal_case1_info.h5', '/node_list');
neighbor_list = h5read('thermal_case1_info.h5', '/neighbor_list');
bound_inner_list = h5read('thermal_case1_info.h5', '/Bound_inner');
bound_outer_list = h5read('thermal_case1_info.h5', '/Bound_outer');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1-2 可视化网格
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1 = figure;
set(gcf, 'unit', 'inch', 'position', [10, 5, 16.00, 8.00]) % figure
patch('Faces', neighbor_list, 'Vertices', node_list, 'FaceColor', 'white');      
hold on
scatter(node_list(:,1), node_list(:,2), 8, 'ko', 'filled')
% figure % 只画边界
% for ii = 1:size(bound_inner_list, 1)
%     index = bound_inner_list(ii, 1);
%     scatter(node_list(index, 1), node_list(index, 2), 8, 'ko', 'filled')
%     hold on
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2 参数设定
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
np_node = size(node_list, 1);
np_unit = size(neighbor_list, 1);

bound2_list = bound_inner_list;
np_bound2 = size(bound_inner_list, 1); % 在内部边界我们施加第二类边界条件
bound3_list = bound_outer_list;
np_bound3 = size(bound3_list, 1); % 在外部边界我们施加第三类边界条件

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3-1 进行有限元计算 ： 矩阵组装
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KG_test = zeros(np_node, np_node);
FG_test = zeros(np_node, 1);
qv = 0; % 内部的热源（ 这里设定为 0 ）
qb = 10; % 边界上的热流密度
beta = 6; % 对流 换热系数
Tc = 25; % 流体温度
Jacobi_mat = zeros(2, 2);
lambda_coeff_mat = diag([1, 1]); % 这里是导热系数，取为1

% 单位节点坐标
unit_node_pos_list = zeros(4, 2);
% 面/线积分所用的点和系数
order_hammer = 1; % Hammer求积的阶数，用于面积分
[hammer_points_list, A_coeff_list] = get_2d_hammer_points(order_hammer);
num_hammer_points = size(hammer_points_list, 1);

% 计算内部单元 形函数 和 形函数偏导 的矩阵 （这里仅考虑平面3节点三角形情形）
std_shape_partial_mat = get_std_shape_partial_mat2(1/3, 1/3, 3);
std_shape_mat = get_std_shape_mat2(1/3, 1/3, 3);

for ii = 1:np_unit
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 得到当前四边形的四个节点的指标和节点上的函数值
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    unit_node_indexes = neighbor_list(ii, :);
    unit_node_pos_list = [node_list(unit_node_indexes(1),:); ...
                          node_list(unit_node_indexes(2),:); ...
                          node_list(unit_node_indexes(3),:);];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3节点三角形积分需要用Hammer积分，而不是Gauss积分
    % 但在这里形函数都是线性函数，所以偏导都是常数，Jacobi矩阵也是常数
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 得到 KS 和 FS
    KS = zeros(3,3);
    FS = zeros(3,1);
    
    % 计算 Hammer point点处的Jacobi矩阵（这里仅考虑平面3节点三角形情形）
    Jacobi_mat = std_shape_partial_mat * unit_node_pos_list;
    det_Jacobi = abs(det(Jacobi_mat));
    J_inv_shape_partial_mat = Jacobi_mat \ std_shape_partial_mat; % J^{-1} * \partial_N
    
    % 累加 KS 和 FS
    KS = KS + A_coeff_list(1) * det_Jacobi * (J_inv_shape_partial_mat)' * lambda_coeff_mat * J_inv_shape_partial_mat;
    FS = FS + A_coeff_list(1) * det_Jacobi * qv * transpose(std_shape_mat);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 将 KS 和 FS 组装进 KG_test 和 FG_test 矩阵中
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for jj = 1:3
        FG_test(unit_node_indexes(jj)) = FG_test(unit_node_indexes(jj)) + FS(jj);
        for kk = 1:3
            KG_test(unit_node_indexes(jj), unit_node_indexes(kk)) = KG_test(unit_node_indexes(jj), unit_node_indexes(kk)) + KS(jj, kk);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3-2 边界条件施加
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 边界节点的组装
% 边界条件1 常值边界条件（这里没有，可以参考test_electric_case1）
% for boundary_index = 1:np_bound1
% end

% 计算边界上的 形函数 和 形函数偏导 的矩阵 （一维 / 高斯点）
order_gauss = 1; % Gauss求积的阶数，用于边界上的线积分
[gauss_points_list, B_coeff_list] = get_gauss_points(order_gauss);
std_1d_shape_mat = [get_std_1d_shape_mat(gauss_points_list(1)), ...
                            get_std_1d_shape_mat(gauss_points_list(2))];

% 边界条件2
% FG_test_bound2 = zeros(np_node, 1);
for boundary_index = 1:np_bound2
    FS_bound2 = zeros(2, 1);
    
    node_left_index = bound2_list(boundary_index, 1);
    node_right_index = bound2_list(boundary_index, 2);
    
    p_left = node_list(node_left_index, :);
    p_right = node_list(node_right_index, :);
    
    % 进行线积分
    % 计算Jacobi (不要忘记除以2)
    det_J = sqrt(sum((p_right - p_left).^2)) / 2; 
    % 由于热源在内部，所以不需要在前面添加负号
    for i = 1:length(gauss_points_list)
        gauss_coeff = B_coeff_list(i);   
        FS_bound2 = FS_bound2 + gauss_coeff * std_1d_shape_mat(:, i);
    end
    FS_bound2 = FS_bound2 * det_J * qb;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 矩阵元累加
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FG_test_bound2(node_left_index) =  FG_test_bound2(node_left_index) + FS_bound2(1);
    % FG_test_bound2(node_right_index) = FG_test_bound2(node_right_index) + FS_bound2(2);
    FG_test(node_left_index) =  FG_test(node_left_index) + FS_bound2(1);
    FG_test(node_right_index) = FG_test(node_right_index) + FS_bound2(2);
end

% 边界条件3
% FG_test_bound3 = zeros(np_node, 1);
% KG_test_bound3 = zeros(np_node, np_node);
for boundary_index = 1:np_bound3
    FS_bound3 = zeros(2, 1);
    KS_bound3 = zeros(2, 2);
    
    node_left_index = bound3_list(boundary_index, 1);
    node_right_index = bound3_list(boundary_index, 2);
    
    p_left = node_list(node_left_index, :);
    p_right = node_list(node_right_index, :);
    
    % 进行线积分
    % 计算Jacobi (不要忘记除以2)
    det_J = sqrt(sum((p_right - p_left).^2)) / 2; 
    for i = 1:length(gauss_points_list)
        gauss_coeff = B_coeff_list(i);
        
        FS_bound3 = FS_bound3 + gauss_coeff * std_1d_shape_mat(:, i);
        KS_bound3 = KS_bound3 + gauss_coeff * std_1d_shape_mat(:, i) * transpose(std_1d_shape_mat(:, i));
    end
    FS_bound3 = FS_bound3 * det_J * beta * Tc;
    KS_bound3 = KS_bound3 * det_J * beta;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 矩阵元累加
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FG_test_bound3(node_left_index) = FG_test_bound3(node_left_index) + FS_bound3(1);
    % FG_test_bound3(node_right_index) = FG_test_bound3(node_right_index) + FS_bound3(2);
    
    % KG_test_bound3(node_left_index, node_left_index) = KG_test_bound3(node_left_index, node_left_index) + KS_bound3(1, 1);
    % KG_test_bound3(node_left_index, node_right_index) = KG_test_bound3(node_left_index, node_right_index) + KS_bound3(1, 2);
    % KG_test_bound3(node_right_index, node_left_index) = KG_test_bound3(node_right_index, node_left_index) + KS_bound3(2, 1);
    % KG_test_bound3(node_right_index, node_right_index) = KG_test_bound3(node_right_index, node_right_index) + KS_bound3(2, 2);
    for ii = 1:2
        index1 = bound3_list(boundary_index, ii);
        FG_test(index1) = FG_test(index1) + FS_bound3(ii);
        for jj = 1:2
            index2 = bound3_list(boundary_index, jj);
            KG_test(index1, index2) = KG_test(index1, index2) + KS_bound3(ii, jj);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3-3 矩阵求解 (求解出来的是温度分布 T(x, y))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U_value_numerical = KG_test \ FG_test;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 后处理 温度分布 T(x, y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 二维伪彩图
fig_value = figure;
patch('Faces', neighbor_list, 'Vertices', node_list, 'FaceVertexCData', U_value_numerical, 'FaceColor', 'interp', 'EdgeAlpha', 1);
colormap(jet)
grid off
shading interp
colorbar;