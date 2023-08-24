function [hammer_points_list, A_coeff_list] = get_2d_hammer_points(order)
    format long  %　调整精度
    switch order
        case 1
            hammer_points_list = [1/3, 1/3];
            A_coeff_list = 1/2;
        case 2
            hammer_points_list = [1/6, 1/6; 2/3, 1/6; 1/6, 2/3];
            A_coeff_list = [1/6; 1/6; 1/6];
        case 3
            hammer_points_list = [1/3, 1/3; 1/5, 1/5; 3/5, 1/5; 1/5, 3/5];
            A_coeff_list = [-27/96; 25/96; 25/96; 25/96];
        otherwise            
            error('Hammer Points only available for m = 1, 2, 3')
    end
end