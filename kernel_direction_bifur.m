function [kernel_matrix_dir_bifur,repeated_input_dir_bifur,target_dir_bifur,theta_initial] = kernel_direction_bifur(theta_initial,eta)

load input_dir_bifur.mat   
load target_dir_bifur.mat

[a,b]=size(input_dir_bifur);
z1(1,:,:) = input_dir_bifur;
repeated_input = repmat(z1, a,1);
% repeated_input =     [input1 input2 input3 ...]
%                       input1 input2 input3 ...
%                       input1 input2 input3 ...
%                         .      .      .

% Parameter Estimation 
kernel_matrix = exp(-theta_initial * sqrt(sum((permute(repeated_input,[2 1 3]) - repeated_input).^2,3)))+.0001;
kernel_matrix_gradient = -(sqrt(sum((permute(repeated_input,[2 1 3]) - repeated_input).^2,3))).*kernel_matrix;
gradient_func(1) = (-1/2)*trace(kernel_matrix\kernel_matrix_gradient) + (1/2)*(target_dir_bifur(:,1)'/kernel_matrix)...
    *kernel_matrix_gradient*(kernel_matrix\target_dir_bifur(:,1));
gradient_func(2) = (-1/2)*trace(kernel_matrix\kernel_matrix_gradient) + (1/2)*(target_dir_bifur(:,2)'/kernel_matrix)...
    *kernel_matrix_gradient*(kernel_matrix\target_dir_bifur(:,2));

% Estimate Theta
while abs(gradient_func(1)) > 5 || abs(gradient_func(2)) > 5
    kernel_matrix = exp(-theta_initial * sqrt(sum((permute(repeated_input,[2 1 3]) - repeated_input).^2,3)))+.0001;
    kernel_matrix_gradient = -(sqrt(sum((permute(repeated_input,[2 1 3]) - repeated_input).^2,3))).*...
        exp(-theta_initial(1) * sqrt(sum((permute(repeated_input,[2 1 3]) - repeated_input).^2,3)));
    if abs(gradient_func(1)) > 5
        gradient_func(1) = (-1/2)*trace(kernel_matrix\kernel_matrix_gradient) + (1/2)*(target_dir_bifur(:,1)'/kernel_matrix)...
            *kernel_matrix_gradient*(kernel_matrix\target_dir_bifur(:,1));
        theta_new = theta_initial + eta*gradient_func(1);
        theta_initial = theta_new;
    end
    
    if abs(gradient_func(2)) > 5
        gradient_func(2) = (-1/2)*trace(kernel_matrix\kernel_matrix_gradient) + (1/2)*(target_dir_bifur(:,2)'/kernel_matrix)...
            *kernel_matrix_gradient*(kernel_matrix\target_dir_bifur(:,2));
        theta_new = theta_initial + eta*gradient_func(2);
        theta_initial = theta_new;
    end
end

kernel_matrix_dir_bifur = exp(-theta_initial * sqrt(sum((permute(repeated_input,[2 1 3]) - repeated_input).^2,3)))+.0001;
repeated_input_dir_bifur = repeated_input(1,:,:);

end