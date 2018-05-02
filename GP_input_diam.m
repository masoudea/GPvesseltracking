function [output] = GP_input_diam(input_image, seed_point_position, mask_size,previous_direction)
%     extract data from image for diameter estimation using radon transform

    zoom = 4;
    [a b] = size(input_image);
    [X,Y] = meshgrid(1:b,1:a);
    [rx,ry] = meshgrid(-mask_size(1):1/zoom:mask_size(1));

    x0 = seed_point_position(1) + ry*cosd(90-previous_direction);
    y0 = seed_point_position(2) + ry*sind(90-previous_direction);

    mask_position_x = x0 + rx*cosd(previous_direction);
    mask_position_y = y0 - rx*sind(previous_direction);

    square_mask = uint8(interp2(X,Y,input_image,mask_position_x,mask_position_y,'cubic'));
    [a,b] = size(square_mask);
    center = (a+1)/2;
    x = 1:a;
    [x,y] = meshgrid(1:a,1:a);
    sig = 10;
    w = exp(-((x-center).^2)/sig);
    weighted_image = double(square_mask).* w;
    radon_in_calculated_direction = sum(weighted_image,2);
    l= length(radon_in_calculated_direction);
    normal_output = interp1(1:l,radon_in_calculated_direction,1:(l-1)/99:l);
    normal_output = normal_output - min(normal_output);
    output = normal_output/max(normal_output);
end