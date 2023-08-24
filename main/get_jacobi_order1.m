%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 计算高斯求积的阶数为 1 时的 Jacobi 矩阵
% 假设 X 方向和 Y 方向取的高斯点数是一样的
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [N_mat_list, Nparial_mat_list] = get_jacobi_order1(gauss_points_list)
    n_gauss_point = length(gauss_points_list);
    Nparial_mat_list = zeros(2, 4, n_gauss_point, n_gauss_point);
    N_mat_list = zeros(4, n_gauss_point, n_gauss_point);
    
    % 计算 每个高斯点处对应的具体的 N函数 的值
    % 计算 每个高斯点处对应的具体的 N函数偏导 的矩阵
    for jj = 1:n_gauss_point
        r_j = gauss_points_list(jj);
        for kk = 1:n_gauss_point
            s_k = gauss_points_list(kk);
            % 计算 每个高斯点处对应的具体的 N函数 的值
            N1 = 0.25 * (1 - r_j) * (1 - s_k);
            N2 = 0.25 * (1 + r_j) * (1 - s_k);
            N3 = 0.25 * (1 + r_j) * (1 + s_k);
            N4 = 0.25 * (1 - r_j) * (1 + s_k);
            N_mat_list(:, jj, kk) = [N1, N2, N3, N4];
            
            % 计算 每个高斯点处对应的具体的 N函数偏导 的矩阵
            Nparial_mat_list(1, 1, jj, kk) = - 0.25 * (1 - s_k);
            Nparial_mat_list(2, 1, jj, kk) = - 0.25 * (1 - r_j);
            Nparial_mat_list(1, 2, jj, kk) =   0.25 * (1 - s_k);
            Nparial_mat_list(2, 2, jj, kk) = - 0.25 * (1 + r_j);
            Nparial_mat_list(1, 3, jj, kk) =   0.25 * (1 + s_k);
            Nparial_mat_list(2, 3, jj, kk) =   0.25 * (1 + r_j);
            Nparial_mat_list(1, 4, jj, kk) = - 0.25 * (1 + s_k);
            Nparial_mat_list(2, 4, jj, kk) =   0.25 * (1 - r_j);
        end
    end

end