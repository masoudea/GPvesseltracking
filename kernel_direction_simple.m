function [kernel_matrix_dir_simple,repeated_input_dir_simple,target_dir_simple,theta_initial] = kernel_direction_simple(theta_initial,eta)

load input_dir_simple.mat   
load target_dir_simple.mat

[a,b]=size(input_dir_simple);
z1(1,:,:) = input_dir_simple;
repeated_input = repmat(z1, a,1);
% repeated_input =     [input1 input2 input3 ...]
%                       input1 input2 input3 ...
%                       input1 input2 input3 ...
%                         .      .      .

% Parameter Estimation 
kernel_matrix = exp(-theta_initial * sqrt(sum((permute(repeated_input,[2 1 3]) - repeated_input).^2,3)))+.0001;
kernel_matrix_gradient = -(sqrt(sum((permute(repeated_input,[2 1 3]) - repeated_input).^2,3))).*kernel_matrix;
gradient_func = (-1/2)*trace(kernel_matrix\kernel_matrix_gradient) + (1/2)*(target_dir_simple'/kernel_matrix)...
    *kernel_matrix_gradient*(kernel_matrix\target_dir_simple);

% Estimate Theta
while abs(gradient_func) > 5
    kernel_matrix = exp(-theta_initial * sqrt(sum((permute(repeated_input,[2 1 3]) - repeated_input).^2,3)))+.0001;
    kernel_matrix_gradient = -(sqrt(sum((permute(repeated_input,[2 1 3]) - repeated_input).^2,3))).*...
        exp(-theta_initial * sqrt(sum((permute(repeated_input,[2 1 3]) - repeated_input).^2,3)));

    gradient_func = (-1/2)*trace(kernel_matrix\kernel_matrix_gradient) + (1/2)*(target_dir_simple(:,1)'/kernel_matrix)...
        *kernel_matrix_gradient*(kernel_matrix\target_dir_simple(:,1));
    theta_new = theta_initial + eta*gradient_func(1);
    theta_initial = theta_new;
end

kernel_matrix_dir_simple = exp(-theta_initial * sqrt(sum((permute(repeated_input,[2 1 3]) - repeated_input).^2,3)))+.0001;
repeated_input_dir_simple = repeated_input(1,:,:);

end