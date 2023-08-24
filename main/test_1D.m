%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test demo in 1D finite element analysis
% (\partial^2 / \partial x^2)u(x) = 2 
% with boundary conditions u(1) = u(-1) = 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 输入参数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%x%%%%%%%%%%%
format long
x_L_end = -1.0;
x_R_end = 1.0;
n_unit = 150;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 单元划分
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
step_length = (x_R_end - x_L_end) / n_unit;
NE = zeros(n_unit, 2);
node_pos_list = zeros(n_unit + 1, 1); % 每个节点的位置坐标
for i = 1:(n_unit + 1)
    if i <= n_unit
        NE(i, 1) = i;
        NE(i, 2) = i + 1;
    end
    node_pos_list(i) = x_L_end + step_length * (i - 1);
end

order_gauss = 1; % 高斯求积的阶数
[gauss_points_list, A_coeff_list] = get_gauss_points(order_gauss);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 矩阵组装 （in the bulk）：非边界点
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KG = zeros(n_unit + 1, n_unit + 1);
FG = zeros(n_unit + 1, 1);
KS = zeros(2, 2);
FS = zeros(2, 1);
% 单位元的形函数 ( shape function )
shape_func_matrix = zeros(order_gauss + 1, order_gauss + 1); % 等参变换
shape_func_partial_matrix = zeros(order_gauss + 1, order_gauss + 1);
% 对每个高斯点进行计算
for ii = 1:(order_gauss + 1)
    h1 = 1 / 2 * (1 - gauss_points_list(ii)); % h1 = 1 / 2 * (1 - r)
    h2 = 1 / 2 * (1 + gauss_points_list(ii)); % h2 = 1 / 2 * (1 + r)
    shape_func_matrix(1, ii) = h1;
    shape_func_matrix(2, ii) = h2;
    shape_func_partial_matrix(1, ii) = - 0.5;
    shape_func_partial_matrix(2, ii) = 0.5;
    
    KS = KS - shape_func_partial_matrix(:, ii) * shape_func_partial_matrix(:, ii)';
    FS = FS + 2 * shape_func_matrix(:, ii);
end

for ii = 1:n_unit
    index_L_end = NE(ii, 1);
    index_R_end = NE(ii, 2);
    
    x_L = node_pos_list(index_L_end);
    x_R = node_pos_list(index_R_end);
    
    % length of the range
    Len_range = (x_R - x_L);
    % Jacobi
    Jacobi = 1 / 2 * (x_R - x_L);
    
    KG(ii:ii+1, ii:ii+1) = KG(ii:ii+1, ii:ii+1) + KS * Len_range / 2 / (Jacobi)^2;
    FG(ii:ii+1) = FG(ii:ii+1) + FS * Len_range / 2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 矩阵组装 （along the edge）: 施加边界条件
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KG(1,1) = 1;
KG(1,2:end) = 0;
KG(end,1:end-1) = 0;
KG(end,end) = 1;

FG(1) = 1;
FG(end) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 求解函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = KG \ FG;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 后处理(作图)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
% 数值解
plot(node_pos_list, U, 'r')
hold on
% 解析解
x_exact = linspace(-1,1,201);
y_exact = x_exact .^ 2;
plot(x_exact, y_exact, 'r--')