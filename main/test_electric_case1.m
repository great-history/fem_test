%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1-1 读取网格(这里采用四边形划分)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
node_list = h5read('electric_case1_info.h5', '/node_list');
neighbor_list = h5read('electric_case1_info.h5', '/neighbor_list');
% 边界节点：第1列存放节点指标，第二列存放节点上的值
bound_inner_list = h5read('electric_case1_info.h5', '/Bound_inner');
bound_outer_list = h5read('electric_case1_info.h5', '/Bound_outer');
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

bound_inner_value = 500;
bound_inner_list(:, 2) = bound_inner_value;
bound_outer_value = 0;
bound_outer_list(:, 2) = bound_outer_value;
boundary_list = [bound_inner_list; bound_outer_list];
np_boundary = size(boundary_list, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3-1 进行有限元计算 ： 矩阵组装
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KG = zeros(np_node, np_node);
FG = zeros(np_node, 1);
QF = 0; % 不考虑电荷分布
Jacobi_mat = zeros(2, 2);
material_coeff_mat = diag([100, 100]); % 这里是介电常数，取为100

% 内部节点的组装
unit_node_pos_list = zeros(4, 2);
% 计算 形函数 和 形函数偏导 的矩阵
order_gauss = 1; % 高斯求积的阶数
[gauss_points_list, A_coeff_list] = get_gauss_points(order_gauss);
[N_mat_list, Nparial_mat_list] = get_jacobi_order1(gauss_points_list);

for ii = 1:np_unit
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 得到当前四边形的四个节点的指标和节点上的函数值
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    unit_node_indexes = neighbor_list(ii, :);
    unit_node_pos_list = [node_list(unit_node_indexes(1),:);node_list(unit_node_indexes(2),:); ...
                            node_list(unit_node_indexes(3),:);node_list(unit_node_indexes(4),:)];
    % 得到 KS 和 FS
    KS = zeros(4,4);
    FS = zeros(4,1);
    for jj = 1:2
        r_j = gauss_points_list(jj);
        A_j = A_coeff_list(jj);
        for kk = 1:2
            s_k = gauss_points_list(kk);
            A_k = A_coeff_list(kk);
            
            Jacobi_mat = squeeze(Nparial_mat_list(:, :, jj, kk)) * unit_node_pos_list; % 该高斯点处的 Jacobi 矩阵
            det_Jacobi_mat = abs(det(Jacobi_mat));
            J_inv_N_partial_mat = Jacobi_mat \ squeeze(Nparial_mat_list(:, :, jj, kk)); % J^{-1} * \partial_N
            
            % 累加 KS 和 FS
            KS = KS - A_j * A_k * det_Jacobi_mat * (J_inv_N_partial_mat)' * material_coeff_mat * J_inv_N_partial_mat;
            FS = FS + A_j * A_k * det_Jacobi_mat * QF * squeeze(N_mat_list(:, jj, kk));
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 将 KS 和 FS 组装进 KG 和 FG 矩阵中
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for jj = 1:4
        FG(unit_node_indexes(jj)) = FG(unit_node_indexes(jj)) + FS(jj);
        for kk = 1:4
            KG(unit_node_indexes(jj), unit_node_indexes(kk)) = KG(unit_node_indexes(jj), unit_node_indexes(kk)) + KS(jj, kk);
        end
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3-2 边界条件施加
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 边界节点的组装
for boundary_index = 1:np_boundary
    node_index = boundary_list(boundary_index, 1);
    node_value = boundary_list(boundary_index, 2);
    KG(node_index, :) = 0;
    KG(node_index, node_index) = 1;
    FG(node_index) = node_value;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3-3 矩阵求解 (求解出来的是电势)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U_value_numerical = KG \ FG;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4 进一步求解电场矢量 (即电势的梯度 U_gradient)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% U_gradient_numerical 第一个维度是节点指标, 
%　第二个维度的维数是3, 分别对应 x方向梯度 / y方向梯度 / 该节点属于几个四边形
U_gradient_numerical = zeros(np_node, 3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4-1 首先求出标准单元的每个节点的形函数偏导
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 节点(-1,-1)
std_shape_partial_mat1 = get_std_shape_partial_mat(-1, -1);
% 节点(1,-1)
std_shape_partial_mat2 = get_std_shape_partial_mat(1, -1);
% 节点(1,1)
std_shape_partial_mat3 = get_std_shape_partial_mat(1, 1);
% 节点(-1,1)
std_shape_partial_mat4 = get_std_shape_partial_mat(-1, 1);

for ii = 1:np_unit % 对每个四边形进行循环
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 得到当前四边形的四个节点的指标和节点上的函数值
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    unit_node_indexes = neighbor_list(ii, :);
    unit_node_pos_list = [node_list(unit_node_indexes(1),:); ...
                          node_list(unit_node_indexes(2),:); ...
                          node_list(unit_node_indexes(3),:); ...
                          node_list(unit_node_indexes(4),:)];
    unit_node_value_list = [U_value_numerical(unit_node_indexes(1)); ...
                            U_value_numerical(unit_node_indexes(2)); ...
                            U_value_numerical(unit_node_indexes(3)); ...
                            U_value_numerical(unit_node_indexes(4));];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 得到当前四边形的四个节点的指标和节点上的函数值
    % 标准四边形单位的四个顶点坐标分别为 (-1,-1) (1,-1) (1,1) (-1,1)
    % 分别计算每个节点对应的Jacobi矩阵和(r,s)空间上的偏导矩阵，然后做矩阵运算
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 节点(-1,-1)
    jacobi_mat = std_shape_partial_mat1 * unit_node_pos_list;
    gradient_rs_mat = std_shape_partial_mat1 * unit_node_value_list;
    gradient_xy_mat = jacobi_mat \ gradient_rs_mat;
    U_gradient_numerical(unit_node_indexes(1), 1:2) = U_gradient_numerical(unit_node_indexes(1), 1:2) + transpose(gradient_xy_mat);
    U_gradient_numerical(unit_node_indexes(1), 3) = U_gradient_numerical(unit_node_indexes(1), 3) + 1;
    
    % 节点(1,-1)
    jacobi_mat = std_shape_partial_mat2 * unit_node_pos_list;
    gradient_rs_mat = std_shape_partial_mat2 * unit_node_value_list;
    gradient_xy_mat = jacobi_mat \ gradient_rs_mat;
    U_gradient_numerical(unit_node_indexes(2), 1:2) = U_gradient_numerical(unit_node_indexes(2), 1:2) + transpose(gradient_xy_mat);
    U_gradient_numerical(unit_node_indexes(2), 3) = U_gradient_numerical(unit_node_indexes(2), 3) + 1;
    
    % 节点(1,1)
    jacobi_mat = std_shape_partial_mat3 * unit_node_pos_list;
    gradient_rs_mat = std_shape_partial_mat3 * unit_node_value_list;
    gradient_xy_mat = jacobi_mat \ gradient_rs_mat;
    U_gradient_numerical(unit_node_indexes(3), 1:2) = U_gradient_numerical(unit_node_indexes(3), 1:2) + transpose(gradient_xy_mat);
    U_gradient_numerical(unit_node_indexes(3), 3) = U_gradient_numerical(unit_node_indexes(3), 3) + 1;
    
    % 节点(-1,1)
    jacobi_mat = std_shape_partial_mat4 * unit_node_pos_list;
    gradient_rs_mat = std_shape_partial_mat4 * unit_node_value_list;
    gradient_xy_mat = jacobi_mat \ gradient_rs_mat;
    U_gradient_numerical(unit_node_indexes(4), 1:2) = U_gradient_numerical(unit_node_indexes(4), 1:2) + transpose(gradient_xy_mat);
    U_gradient_numerical(unit_node_indexes(4), 3) = U_gradient_numerical(unit_node_indexes(4), 3) + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4-2 进行平均化
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:np_node
    num_cum = U_gradient_numerical(ii, 3);
    U_gradient_numerical(ii, 1) = U_gradient_numerical(ii, 1) / num_cum;
    U_gradient_numerical(ii, 2) = U_gradient_numerical(ii, 2) / num_cum;
end

U_gradient_abs_numerical_list = ((U_gradient_numerical(:,1)).^2 + (U_gradient_numerical(:,2)).^2).^0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 后处理 电势 / 电场矢量
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 二维伪彩图
fig_value = figure;
patch('Faces', neighbor_list, 'Vertices', node_list, 'FaceVertexCData', U_value_numerical, 'FaceColor', 'interp', 'EdgeAlpha', 1);
colormap(jet)
grid off
shading interp
colorbar;

% 二维伪彩图
fig_gradient_x = figure;
patch('Faces', neighbor_list, 'Vertices', node_list, 'FaceVertexCData', U_gradient_numerical(:,1), 'FaceColor', 'interp', 'EdgeAlpha', 1);
colormap(jet)
grid off
shading interp
colorbar;

% 二维伪彩图
fig_gradient_y = figure;
patch('Faces', neighbor_list, 'Vertices', node_list, 'FaceVertexCData', U_gradient_numerical(:,2), 'FaceColor', 'interp', 'EdgeAlpha', 1);
colormap(jet)
grid off
shading interp
colorbar;

% 二维伪彩图
fig_gradient_abs = figure;
patch('Faces', neighbor_list, 'Vertices', node_list, 'FaceVertexCData', U_gradient_abs_numerical_list, 'FaceColor', 'interp', 'EdgeAlpha', 1);
colormap(jet)
grid off
shading interp
colorbar;