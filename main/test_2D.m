%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test demo in 2D finite element analysis
% 二维泊松方程: \Delta u = 4, -1 <= x,y <= 1
% 边界条件 : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 输入参数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%x%%%%%%%%%%%
format long
X_L = -1.0;
X_R = 1.0;
Y_U = 1.0;
Y_D = -1.0;
np_x = 21;
np_y = 21;
np_node = np_x * np_y;
np_boundary = 2 * (np_x + np_y) - 4;
np_unit = (np_x - 1) * (np_y - 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 单元划分
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_value_list = linspace(X_L, X_R, np_x);
y_value_list = linspace(Y_D, Y_U, np_y);
node_pos_list = zeros(np_node, 2);
boundary_list = zeros(np_boundary, 2); % 存放边界点的指标和边界值
unit_index_list = zeros(np_unit, 4); % 每个区域四个点的指标

% 得到 node_pos_list / boundary_list
node_index = 0;
boundary_index = 0;
for j = 1:np_y
    for i = 1:np_x
        node_index = node_index + 1;
        x_value = x_value_list(i);
        y_value = y_value_list(j);
        node_pos_list(node_index, :) = [x_value, y_value];
        % 确定边界节点
        if (x_value == X_L) || (x_value == X_R) || (y_value == Y_U) || (y_value == Y_D)
            boundary_index = boundary_index + 1;
            boundary_list(boundary_index, 1) = node_index;
            boundary_list(boundary_index, 2) = x_value^2 + y_value^2;
        end
    end
end

% 得到 unit_index_list
unit_index = 0;
for j = 1:(np_y - 1)
    left_down_index = (j - 1) * np_x;
    for i = 1:(np_x - 1)
        left_down_index = left_down_index + 1;
        unit_index = unit_index + 1;
        % left_down
        unit_index_list(unit_index, 1) = left_down_index;
        % right_down
        unit_index_list(unit_index, 2) = left_down_index + 1;
        % right_up
        unit_index_list(unit_index, 3) = left_down_index + 1 + np_x;
        % left_up
        unit_index_list(unit_index, 4) = left_down_index + np_x;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 画出网格
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
% 撒点
scatter(node_pos_list(:,1),node_pos_list(:,2),'ko','filled');
hold on
% 给每个四边形单元画线
for i=1:np_unit
    p = zeros(5,2);
    for j=1:4
        node_index = unit_index_list(i,j);
        p(j,:) = node_pos_list(node_index, :);
    end
    p(5,:) = node_pos_list(unit_index_list(i,1), :);
    plot(p(:,1),p(:,2),'k-');
    hold on
end

order_gauss = 1; % 高斯求积的阶数
[gauss_points_list, A_coeff_list] = get_gauss_points(order_gauss);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 矩阵组装
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KG = zeros(np_node, np_node);
FG = zeros(np_node, 1);
QF = 4; % 微分方程的源，这里取为常数，也可以是x,y的函数
Jacobi_mat = zeros(2, 2);
material_coeff_mat = diag([1,1]); % 可以是介电常数或者其它物理量
[N_mat_list, Nparial_mat_list] = get_jacobi_order1(gauss_points_list); % 计算 N函数 和 N函数偏导 的矩阵
unit_node_pos_list = zeros(4, 2);

% 内部节点
for ii = 1:np_unit
    unit_node_indexes = unit_index_list(ii, :);
    % 单元上四个节点的坐标
    unit_node_pos_list = [node_pos_list(unit_node_indexes(1),:);node_pos_list(unit_node_indexes(2),:); ...
                            node_pos_list(unit_node_indexes(3),:);node_pos_list(unit_node_indexes(4),:)];
    % 得到 KS 和 FS
    KS = zeros(4,4);
    FS = zeros(4,1);
    for jj = 1:2
        r_j = gauss_points_list(jj);
        A_j = A_coeff_list(jj);
        for kk = 1:2
            s_k = gauss_points_list(kk);
            A_k = A_coeff_list(kk);
            
            Jacobi_mat = Nparial_mat_list(:, :, jj, kk) * unit_node_pos_list; % 该高斯点处的 Jacobi 矩阵
            det_Jacobi_mat = abs(det(Jacobi_mat));
            J_inv_N_partial_mat = Jacobi_mat \ Nparial_mat_list(:, :, jj, kk);
            
            % 累加 KS 和 FS
            KS = KS - A_j * A_k * det_Jacobi_mat * (J_inv_N_partial_mat)' * material_coeff_mat * J_inv_N_partial_mat;
            FS = FS + A_j * A_k * det_Jacobi_mat * QF * N_mat_list(:, jj, kk);
        end
    end
    
    % 将 KS 和 FS 组装进 KG 和 FG 矩阵中
    for jj = 1:4
        FG(unit_node_indexes(jj)) = FG(unit_node_indexes(jj)) + FS(jj);
        for kk = 1:4
            KG(unit_node_indexes(jj), unit_node_indexes(kk)) = KG(unit_node_indexes(jj), unit_node_indexes(kk)) + KS(jj, kk);
        end
    end
    
end

% 边界节点
for boundary_index = 1:np_boundary
    node_index = boundary_list(boundary_index, 1);
    node_value = boundary_list(boundary_index, 2);
    KG(node_index, :) = 0;
    KG(node_index, node_index) = 1;
    FG(node_index) = node_value;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 矩阵求解
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U_numerical = KG \ FG;
U_exact = zeros(np_x, np_y);
for ii = 1:np_x
    x_value = x_value_list(ii);
    for jj = 1:np_y
        y_value = y_value_list(jj);
        U_exact(ii, jj) = x_value^2 + y_value^2;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 后处理
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 数值解
% subplot(2,3,1),POST_FIGURE(unit_index_list, node_pos_list, U_numerical);
% title('2D数值解')
% % subplot(2,3,4),POST_3D_FIGURE(X_L, X_R, Y_D, Y_U, np_x, U_numerical, node_pos_list)
% % title('3D数值解')
% % 准确解
% subplot(2,3,2),POST_FIGURE(unit_index_list, node_pos_list, U_exact);
% title('2D准确解')
% subplot(2,3,5),POST_3D_FIGURE(X_L, X_R, Y_D, Y_U, np_x, U_exact, node_pos_list)
% title('3D准确解')
% 误差
% error=((UZ-U).^2).^0.5;
% subplot(2,3,3),POST_FIGURE(NE,NODE,error);
% title('2D误差分布')
% subplot(2,3,6),POST_3D_FIGURE(X1,step,X2,Y1,Y2,error,NODE)
% title('3D数值解')