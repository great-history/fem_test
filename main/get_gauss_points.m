%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 根据 order 给出对应的高斯点 gauss_points_list 和 系数 A_coeff_list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gauss_points_list, A_coeff_list] = get_gauss_points(order)
    format long  %　调整精度
    switch order
        case 0
            gauss_points_list = 0;
            A_coeff_list = 2;
        case 1
            gauss_points_list = [-1/sqrt(3), 1/sqrt(3)];
            A_coeff_list = [1, 1];
        case 2
            gauss_points_list = [-sqrt(0.6), 0, sqrt(0.6)];
            A_coeff_list = [5/9, 8/9, 5/9];
        case 3
            gauss_points_list = [- sqrt(525 + 70 * sqrt(30)) / 35, - sqrt(525 - 70 * sqrt(30)) / 35, sqrt(525 - 70 * sqrt(30)) / 35, sqrt(525 + 70 * sqrt(30)) / 35];
            A_coeff_list = [(18 - sqrt(30)) / 36, (18 + sqrt(30)) / 36, (18 + sqrt(30)) / 36, (18 - sqrt(30)) / 36];
        case 4
            gauss_points_list = [- sqrt(245 + 14 * sqrt(70)) / 21, - sqrt(245 - 14 * sqrt(70)) / 21, 0, sqrt(245 - 14 * sqrt(70)) / 21, sqrt(245 + 14 * sqrt(70)) / 21];
            A_coeff_list = [(322 - 13*sqrt(70))/900, (322 + 13*sqrt(70))/900, 0, (322 + 13*sqrt(70))/900, (322 - 13*sqrt(70))/900];
        otherwise
            error('Gauss Points only available for m = 0,1,2,3,4')
    end
end