% Tracking and diameter estimation of 2D blood vessels using Gaussian Process 
% 
% This script track center points and diameter of blood vessels, which is 
% an ongoing challenge in medical image analysis. We hypothesize that the 
% curvature and the diameter of blood vessels are Gaussian processes (GPs).
% Local Radon transform, which is robust against noise, is subsequently 
% used to compute the features and train the GPs. By learning the 
% kernelized covariance matrix from training data, vessel direction and 
% its diameter are estimated. In order to detect bifurcations, multiple 
% GPs are used and the difference between their corresponding predicted 
% directions is quantified. 
% 
% References: 
% Masoud Elhami Asl, et al. "Tracking and diameter estimation of retinal 
% vessels using Gaussian process and Radon transform." Journal of Medical 
% Imaging 4.3 (2017): 034006.
% 
% This algorithm is the result of many hours of work and problem solving.
% Please cite the above paper in case you find the script useful in your 
% own research. I might be able to help you with your queries, if you send 
% me (Masoud Elhami Asl) an email to: masoudelhamiasl@gmail.com
% 
% Developed and Copyrighted by Masoud Elhami Asl (2017)
%

clear all
close all
clc

%% parameters initialization
step_size = 1;
mask_size_dir = 12; 
mask_size_diam = 6;
theta_dir_simple = 0.2655; % kernel matrix
theta_dir_bifur = 0.3395; % kernel matrix
theta_diameter = 0.07; % kernel_diameter matrix
eta1 = 1e-6 ; % parameter estimation and gradient desent
eta2 = 1e-6; % parameter estimation and gradient desent
eta_diameter = 1e-5;

centerpoints = zeros(1,2);
vessel_edge = zeros(1,4);
bifur = 0;
diameter = 0;

%% channel G in RGB color space is selected
Image = imread('16_test.tif');
input_image = double(Image(:,:,2));

%% seed point selection
centerpoints(1,1) = 455;
centerpoints(1,2) = 317;
direction = -125;
diameter = 5;

% uncomment this part to select desired seed point manually
% disp('select two initial points inside the vessel you want to track and')
% disp('two points on the edge of blood vessel to compute initial diameter')
% disp('press Enter to continue')
% pause
% [seedpoint_position, initial_direction, vessel_initial_diameter] = seed_point_selection(input_image);
% centerpoints(1,1) = seedpoint_position(1,1);
% centerpoints(1,2) = seedpoint_position(1,2);
% direction = initial_direction;
% diameter = vessel_initial_diameter;

vessel_edge(end,1)= centerpoints(1,1) + (diameter/2) * cosd(direction+90);
vessel_edge(end,2)= centerpoints(1,2) - (diameter/2) * sind(direction+90);
vessel_edge(end,3)= centerpoints(1,1) + (diameter/2) * cosd(direction-90);
vessel_edge(end,4)= centerpoints(1,2) - (diameter/2) * sind(direction-90);

%% load training data and compute a 3D matrix which is used to generate kernel matrix without any 'for' loop
disp('GP Hyperparameter Estimating')
disp('1 out of 3')
[kernel_matrix_dir_simple,repeated_input_dir_simple,target_dir_simple,theta_dir_simple]= kernel_direction_simple(theta_dir_simple,eta1);
disp('2 out of 3')
[kernel_matrix_dir_bifur,repeated_input_dir_bifur,target_dir_bifur,theta_dir_bifur]= kernel_direction_bifur(theta_dir_bifur,eta2);
disp('3 out of 3')
[kernel_matrix_diameter,repeated_input_diameter,target_diameter,theta_diameter] = kernel_diameter(theta_diameter,eta_diameter);

disp('Tracking Process Started')
for j=1:400
    
    if mod(j,50) == 0
        fprintf('%d out of 415 \n', j)
    end
    
    %% simple vessel tracking  
    new_input = GP_input_dir(input_image,centerpoints(end,:),mask_size_dir,direction);

    [a,b,c] = size(repeated_input_dir_simple);
    z1(1,1,:) = new_input;
    repeated_new_input = repmat(z1,1,b);
    diff_new_input_and_other_input = (repeated_input_dir_simple - repeated_new_input).^2;
    norm_diff_new_other(:,:,1) = -theta_dir_simple * sqrt(sum(diff_new_input_and_other_input,3));

    kernel_kxsx = exp(sum(norm_diff_new_other,3))';
    kernel_kxsxs= 1.0001;
    mx_xnew_dir_simple= kernel_kxsx' / kernel_matrix_dir_simple * target_dir_simple;
    direction = direction + mx_xnew_dir_simple;
    direction = check_direction(direction);

    %% bifurcation detection         
    [a,b,c] = size(repeated_input_dir_bifur);
    repeated_new_input = repmat(z1,1,b);

    diff_new_input_and_other_input = (repeated_input_dir_bifur - repeated_new_input).^2;
    norm_diff_new_other2(:,:,1) = -theta_dir_bifur * sqrt(sum(diff_new_input_and_other_input,3));

    kernel_kxsx = exp(sum(norm_diff_new_other2,3))';
    kernel_kxsxs= 1.0001;
    mx_xnew_bifur(1)= kernel_kxsx' / kernel_matrix_dir_bifur * target_dir_bifur(:,1);
    mx_xnew_bifur(2)= kernel_kxsx' / kernel_matrix_dir_bifur * target_dir_bifur(:,2);
    bifur(end+1) = abs(mx_xnew_bifur(1) - mx_xnew_bifur(2));

    %% diameter Estimation 
    new_input_diameter = GP_input_diam(input_image,centerpoints(end,:),mask_size_diam,direction);
    [a,b,c] = size(repeated_input_diameter);
    z2(1,1,:) = new_input_diameter;
    repeated_new_input_diameter = repmat(z2,1,b);

    diff_new_input_diameter_and_other = (repeated_input_diameter - repeated_new_input_diameter).^2;
    norm_diff_new_input_diameter_with_other = -theta_diameter*sqrt(sum(diff_new_input_diameter_and_other,3));
    kernel_kxsx_diameter = exp(sum(norm_diff_new_input_diameter_with_other,3))';
    kernel_kxsxs_diameter= 1.0001;
    mx_xnew_diameter = diameter + kernel_kxsx_diameter' /kernel_matrix_diameter * (target_diameter - diameter);
    diameter = max(kernel_kxsx) * mx_xnew_diameter;

    centerpoints(end+1,1) = centerpoints(end,1) + (step_size*cosd(direction));
    centerpoints(end,2) = centerpoints(end-1,2) - (step_size*sind(direction));

    vessel_edge(end+1,1) = centerpoints(end,1) + (diameter/2) * cosd(direction+90);
    vessel_edge(end,2) = centerpoints(end,2) - (diameter/2) * sind(direction+90);
    vessel_edge(end,3) = centerpoints(end,1) + (diameter/2) * cosd(direction-90);
    vessel_edge(end,4) = centerpoints(end,2) - (diameter/2) * sind(direction-90);    
end

figure(1); imshow(uint8(Image))
hold on
plot(centerpoints(1:4:end,1),centerpoints(1:4:end,2),'w.')
plot(vessel_edge(1:4:end,1),vessel_edge(1:4:end,2),'b.')
plot(vessel_edge(1:4:end,3),vessel_edge(1:4:end,4),'b.')

bf = (bifur(2:end-1) - bifur(1:end-2)) .* (bifur(2:end-1) - bifur(3:end));
bfx = find((bf > 0 & bifur(2:end-1)> 50))+ 1;
figure(2);
p1 = plot(bifur);
hold on
p2 = plot(bfx, bifur(bfx),'r*');
xlabel('center point index')
ylabel('difference between m1 and m2')
legend([p2],{'bifurcation'})
