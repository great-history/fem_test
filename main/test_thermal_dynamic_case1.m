%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 读取网格(这里采用三角形划分)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
node_list = h5read('thermal_dynamic_case1_info.h5', '/node_list');
neighbor_list = h5read('thermal_dynamic_case1_info.h5', '/neighbor_list');
bound_source_list = h5read('thermal_dynamic_case1_info.h5', '/Bound_source');
bound_cond_list = h5read('thermal_dynamic_case1_info.h5', '/Bound_cond');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1-2 可视化网格
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1 = figure;
set(gcf, 'unit', 'inch', 'position', [10, 5, 16.00, 8.00]) % figure
patch('Faces', neighbor_list, 'Vertices', node_list, 'FaceColor', 'white');      
hold on
scatter(node_list(:,1), node_list(:,2), 8, 'ko', 'filled')
% 
% figure % 只画热源的边界
% for ii = 1:size(bound_source_list, 1)
%     index = bound_source_list(ii, 1);
%     scatter(node_list(index, 1), node_list(index, 2), 8, 'ko', 'filled')
%     hold on
% end
% 
% figure % 只画热传导的边界
% for ii = 1:size(bound_cond_list, 1)
%     index = bound_cond_list(ii, 1);
%     scatter(node_list(index, 1), node_list(index, 2), 8, 'ko', 'filled')
%     hold on
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2-1 基本参数设定
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
np_node = size(node_list, 1);
np_unit = size(neighbor_list, 1);

np_bound_source = size(bound_source_list, 1); % 在内部边界我们施加第二类边界条件
np_bound_cond = size(bound_cond_list, 1); % 在外部边界我们施加第三类边界条件

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2-2 微分方程中的相关系数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
density = 7.85E-9; % 材料的密度，单位是 tonne / mm^3
heat_capacity = 4.36E8; % 材料热容，单位是 mJ / (tonne * K)
lambda_x = 60.5; % 材料沿x方向热传导系数，单位是 W / (m * C)
lambda_y = 60.5; % 材料沿y向热传导系数，单位是 W / (m * C)
beta = 1; % 与外界的热交换系数 6
temp_ext = 25; % 外界温度
power_source = 100; % 热源功率（ 这里设定为 100 W ）
power_material = 0; % 材料自身的发热功率 （这里设为 0 W）

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3-1 计算所需要用到的矩阵
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K_cond = zeros(np_node, np_node); % C for conductive
F_cond = zeros(np_node, 1);
K_capacity = zeros(np_node, np_node); % P for power

% 计算 形函数 和 形函数偏导 的矩阵
order_gauss = 1; % 高斯求积的阶数
[gauss_points_list, A_coeff_list] = get_gauss_points(order_gauss);
[std_shape_mat, std_shape_partial_mat] = get_jacobi_order1(gauss_points_list);

% Jacobi矩阵
Jacobi_mat = zeros(2, 2);
lambda_coeff_mat = diag([lambda_x, lambda_y]); % 这里是导热系数

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3-2 体内的有限元矩阵
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:np_unit
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 得到当前四边形的四个节点的指标和节点上的函数值
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    unit_node_indexes = neighbor_list(ii, :);
    unit_node_pos_list = [node_list(unit_node_indexes(1),:); ...
                          node_list(unit_node_indexes(2),:); ...
                          node_list(unit_node_indexes(3),:); ...
                          node_list(unit_node_indexes(4),:)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 计算 K_cond / F_cond 的矩阵
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 得到 K_unit 和 FS
    K_unit = zeros(4, 4);
    F_unit = zeros(4, 1);
    
    % 累加 K_unit 和 F_unit
    for jj = 1:2
        r_j = gauss_points_list(jj);
        A_j = A_coeff_list(jj);
        for kk = 1:2
            s_k = gauss_points_list(kk);
            A_k = A_coeff_list(kk);
            
            % 计算 Gauss point点处的Jacobi矩阵 及 其行列式的绝对值
            Jacobi_mat = squeeze(std_shape_partial_mat(:, :, jj, kk)) * unit_node_pos_list; % 该高斯点处的 Jacobi 矩阵
            det_Jacobi_mat = abs(det(Jacobi_mat));
            J_inv_N_partial_mat = Jacobi_mat \ squeeze(std_shape_partial_mat(:, :, jj, kk)); % J^{-1} * \partial_N
            
            % 累加 K_unit 和 F_unit
            K_unit = K_unit + A_j * A_k * det_Jacobi_mat * (J_inv_N_partial_mat)' * lambda_coeff_mat * J_inv_N_partial_mat;
            F_unit = F_unit + A_j * A_k * det_Jacobi_mat * power_material * squeeze(std_shape_mat(:, jj, kk));
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 将 K_unit 和 F_unit 组装进 K_cond 和 F_cond 矩阵中
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for jj = 1:4
        F_cond(unit_node_indexes(jj)) = F_cond(unit_node_indexes(jj)) + F_unit(jj);
        for kk = 1:4
            K_cond(unit_node_indexes(jj), unit_node_indexes(kk)) = K_cond(unit_node_indexes(jj), unit_node_indexes(kk)) + K_unit(jj, kk);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 计算 K_capacity / FP 的矩阵
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 得到 K_unit 和 F_unit
    K_unit = zeros(4, 4);
    F_unit = zeros(4, 1);
    
    % 累加 K_unit 和 F_unit
    for jj = 1:2
        r_j = gauss_points_list(jj);
        A_j = A_coeff_list(jj);
        for kk = 1:2
            s_k = gauss_points_list(kk);
            A_k = A_coeff_list(kk);
            
            % 计算 Gauss point点处的Jacobi矩阵 及 其行列式的绝对值
            Jacobi_mat = squeeze(std_shape_partial_mat(:, :, jj, kk)) * unit_node_pos_list; % 该高斯点处的 Jacobi 矩阵
            det_Jacobi_mat = abs(det(Jacobi_mat));
            % J_inv_N_partial_mat = Jacobi_mat \ squeeze(std_shape_partial_mat(:, :, jj, kk)); % J^{-1} * \partial_N
            
            % 累加 K_unit 和 F_unit
            K_unit = K_unit + A_j * A_k * det_Jacobi_mat * std_shape_mat(:, jj, kk) * transpose(std_shape_mat(:, jj, kk));
        end
    end
    K_unit = K_unit * density * heat_capacity;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 将 K_unit 和 F_unit 组装进 K_capacity 和 FP 矩阵中
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for jj = 1:4
        for kk = 1:4
            K_capacity(unit_node_indexes(jj), unit_node_indexes(kk)) = K_capacity(unit_node_indexes(jj), unit_node_indexes(kk)) + K_unit(jj, kk);
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3-3 边界上的有限元矩阵
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 边界节点的组装
% 边界条件1 常值边界条件（这里没有，可以参考test_electric_case1）

% 计算边界上的 形函数 和 形函数偏导 的矩阵 （一维 / 高斯点）
order_gauss = 1; % Gauss求积的阶数，用于边界上的线积分
[gauss_points_list, B_coeff_list] = get_gauss_points(order_gauss);
std_1d_shape_mat = [get_std_1d_shape_mat(gauss_points_list(1)), ...
                    get_std_1d_shape_mat(gauss_points_list(2))];

% 边界条件2
F_bd_src = zeros(np_node, 1);
for boundary_index = 1:np_bound_source
    F_unit_bound2 = zeros(2, 1);
    
    node_left_index = bound_source_list(boundary_index, 1);
    node_right_index = bound_source_list(boundary_index, 2);
    
    p_left = node_list(node_left_index, :);
    p_right = node_list(node_right_index, :);
    
    % 进行线积分
    % 计算Jacobi (不要忘记除以2)
    det_J = sqrt(sum((p_right - p_left).^2)) / 2; 
    % 由于热源在内部，所以不需要在前面添加负号
    for i = 1:length(gauss_points_list)
        gauss_coeff = B_coeff_list(i);   
        F_unit_bound2 = F_unit_bound2 + gauss_coeff * std_1d_shape_mat(:, i);
    end
    F_unit_bound2 = F_unit_bound2 * det_J * power_source;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 矩阵元累加
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FG_test_bound2(node_left_index) =  FG_test_bound2(node_left_index) + F_unit_bound2(1);
    % FG_test_bound2(node_right_index) = FG_test_bound2(node_right_index) + F_unit_bound2(2);
    F_bd_src(node_left_index) =  F_bd_src(node_left_index) + F_unit_bound2(1);
    F_bd_src(node_right_index) = F_bd_src(node_right_index) + F_unit_bound2(2);
end                

% 边界条件3
F_bd_cond = zeros(np_node, 1);
K_bd_cond = zeros(np_node, np_node);
for boundary_index = 1:np_bound_cond
    F_unit_bound3 = zeros(2, 1);
    K_unit_bound3 = zeros(2, 2);
    
    node_left_index = bound_cond_list(boundary_index, 1);
    node_right_index = bound_cond_list(boundary_index, 2);
    
    p_left = node_list(node_left_index, :);
    p_right = node_list(node_right_index, :);
    
    % 进行线积分
    % 计算Jacobi (不要忘记除以2)
    det_J = sqrt(sum((p_right - p_left).^2)) / 2; 
    for i = 1:length(gauss_points_list)
        gauss_coeff = B_coeff_list(i);
        
        F_unit_bound3 = F_unit_bound3 + gauss_coeff * std_1d_shape_mat(:, i);
        K_unit_bound3 = K_unit_bound3 + gauss_coeff * std_1d_shape_mat(:, i) * transpose(std_1d_shape_mat(:, i));
    end
    F_unit_bound3 = F_unit_bound3 * det_J * beta * temp_ext;
    K_unit_bound3 = K_unit_bound3 * det_J * beta;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 矩阵元累加
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ii = 1:2
        index1 = bound_cond_list(boundary_index, ii);
        F_bd_cond(index1) = F_bd_cond(index1) + F_unit_bound3(ii);
        for jj = 1:2
            index2 = bound_cond_list(boundary_index, jj);
            K_bd_cond(index1, index2) = K_bd_cond(index1, index2) + K_unit_bound3(ii, jj);
        end
    end
end

KG1 = K_bd_cond + K_cond;
FG1 = F_bd_src + F_bd_cond; % F_bd_src 不用加负号

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4 ： 按照时间步进行迭代计算
% 已经验证 ：当ttol = 0.000001时，100步内就能收敛
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4-1 设置收敛参数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ttol = 0.0000001; % 收敛值，当温度变化比例小于此时停止循环
time_step = 1; % 初始时间步，上一时刻与下一时刻之间的步长
time_now = 0; % 当前的时刻

count = 0;
ttol_now = 100;
num_print = 10; % 每隔 num_print 打印一下当前的状态

K_mat = zeros(np_node, np_node);
F_vec = zeros(np_node, 1);
T_last_vec = zeros(np_node, 1); % 初始温度分布，这里全设为0

h5file_name = ['thermal_dynamic_case1_', datestr(now,'yyyy-mm-dd_HH-MM-SS' ), '.h5'];
while ttol_now >= ttol
    count = count + 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 有限元计算
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    K_mat = KG1 * time_step + K_capacity;
    F_vec = FG1 * time_step + K_capacity * T_last_vec;
    T_now_vec = K_mat \ F_vec;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 收敛程度 及 时间步更新
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ttol_now = max(abs((T_now_vec - T_last_vec) ./ T_last_vec));
    T_last_vec = T_now_vec;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 打印信息并保存
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~mod(count, num_print)
        % disp("当前运行到第？步，收敛程度为？")
        fprintf('当前运行到第%d步，收敛程度为%f\n', count, ttol_now);
        h5create(h5file_name, ['/T_vec_', num2str(count)], np_node, 'Datatype', 'double');
        h5write(h5file_name, ['/T_vec_', num2str(count)], T_last_vec);
    end
    
end

if ~mod(count, num_print)    
    num_figs = floor(count / num_print);
else
    fprintf('当前运行到第%d步，收敛程度为%f\n', count, ttol_now);
    h5create(h5file_name, ['/T_vec_', num2str(count)], np_node, 'Datatype', 'double');
    h5write(h5file_name, ['/T_vec_', num2str(count)], T_last_vec);
    
    num_figs = floor(count / num_print) + 1;
end

% h5disp(h5file_name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 后处理 温度分布 T(x, y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 二维伪彩图
fig_value = figure;
patch('Faces', neighbor_list, 'Vertices', node_list, 'FaceVertexCData', T_now_vec, 'FaceColor', 'interp', 'EdgeAlpha', 1);
colormap(jet)
grid off
shading interp
colorbar;

% 二维动态图
% for i=1:num_figs
%     T_now = h5read(h5file_name, ['/T_vec_', num2str(i * num_print)]);
%     
%     patch('Faces', neighbor_list, 'Vertices', node_list, 'FaceVertexCData', T_now, 'FaceColor', 'interp', 'EdgeAlpha', 0);
%     colorbar
%     colormap("jet")
%     axis equal
%     if i~=1
%         a=min(T_now);
%         b=max(T_now);
%         caxis([a,b])
%     end
%     drawnow
%     pause(0.5)
% end