function [kernel_matrix_diameter,repeated_input_diameter,target_diameter,theta_initial] = kernel_diameter(theta_initial,eta_diameter)

load input_diameter.mat   
load target_diameter.mat

[a,b]=size(input_diameter);
z2(1,:,:) = input_diameter;
repeated_input_diameter = repmat(z2, a,1);
% repeated_input_diameter =     [input1 input2 input3 ...]
%                                input1 input2 input3 ...
%                                input1 input2 input3 ...
%                                   .      .      .

% Parameter Estimation for Diameter
kernel_matrix_diameter = exp(-theta_initial(1) * sqrt(sum((permute(repeated_input_diameter,[2 1 3]) - repeated_input_diameter).^2,3)))+.0001;
kernel_matrix_diameter_gradient = -sqrt(sum((permute(repeated_input_diameter,[2 1 3]) - repeated_input_diameter).^2,3)).*kernel_matrix_diameter;
gradient_diameter_func = (-1/2)*trace(kernel_matrix_diameter\kernel_matrix_diameter_gradient) + (1/2)*(target_diameter'/kernel_matrix_diameter)...
    *kernel_matrix_diameter_gradient*(kernel_matrix_diameter\target_diameter);

while abs(gradient_diameter_func) > 5
    kernel_matrix_diameter = exp(-theta_initial(1) * sqrt(sum((permute(repeated_input_diameter,[2 1 3]) - repeated_input_diameter).^2,3)))+.0001;
    kernel_matrix_diameter_gradient = -sqrt(sum((permute(repeated_input_diameter,[2 1 3]) - repeated_input_diameter).^2,3)).*kernel_matrix_diameter;
    gradient_diameter_func = (-1/2)*trace(kernel_matrix_diameter\kernel_matrix_diameter_gradient) + (1/2)*(target_diameter'/kernel_matrix_diameter)...
        *kernel_matrix_diameter_gradient*(kernel_matrix_diameter\target_diameter);
    if abs(gradient_diameter_func) > 5
        theta_new_diameter = theta_initial(1) + eta_diameter*gradient_diameter_func;
        theta_initial(1) = theta_new_diameter;
    end
end

kernel_matrix_diameter = exp(-theta_initial * sqrt(sum((permute(repeated_input_diameter,[2 1 3]) - repeated_input_diameter).^2,3)))+.0001;
repeated_input_diameter = repeated_input_diameter(1,:,:);
end