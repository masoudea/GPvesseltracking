function [output] = GP_input_dir(input_image, seed_point_position, vessel_radius,previous_direction)
%     extract data from image for vessel tracking using radon transform    

    zoom = 4;
    [a b] = size(input_image);
    [X,Y] = meshgrid(1:b,1:a);

    i=1;
    [grid_xx,grid_yy] = meshgrid(seed_point_position(1)-vessel_radius(i):1/zoom:seed_point_position(1)+vessel_radius(i) , ...
        seed_point_position(2)-vessel_radius(i):1/zoom:seed_point_position(2)+vessel_radius(i));
    square_mask = uint8(interp2(X,Y,input_image,grid_xx,grid_yy,'cubic'));
    square_mask = double(square_mask);

    min_mask = min(square_mask(:))/255;
    max_mask = max(square_mask(:))/255;
    square_mask = double(imadjust(uint8(square_mask),[min_mask max_mask], [0 1]));
    Mask_Size = size(square_mask);
    Mask_center = (Mask_Size(1)+1)/2;
    [yy,xx] = ndgrid( (1:Mask_Size(1)),(1:Mask_Size(2)));
    circle_mask= double((xx-Mask_center).^2+(yy-Mask_center).^2 <= (zoom*vessel_radius(i))^2);
    image_circle_mask = square_mask.*circle_mask;
    sigma_2= -Mask_center^2/log(.15);
    [x,y] = meshgrid(-Mask_center+1:1:Mask_center-1);
    z= exp(-((x).^2)/sigma_2 - ((y).^2)/sigma_2);
    circle_mask_gaussian= z.*image_circle_mask;
    rotated_circle_mask = imrotate(circle_mask_gaussian,-previous_direction);

    [a b]= size(rotated_circle_mask);
    seed_point_center = (a+1)/2;
    mask_radius = seed_point_center - 1;
    [XX,YY] = meshgrid(-mask_radius:mask_radius);

    integral_radius = [0: mask_radius/100:mask_radius];
    theta_int = [89:-1:-89];
    grid_x(theta_int+90,:) = cosd(theta_int)'*integral_radius;
    grid_y(theta_int+90,:) = sind(theta_int)'*integral_radius;

    ZI = interp2(XX,YY,rotated_circle_mask,grid_x,grid_y,'cubic');
    integral_179 = sum(ZI');

    integral_179 = integral_179 - min(integral_179);
    integral_179_normal = integral_179/max(integral_179);
    output((i-1)*179+1:i*179) = integral_179(179:-1:1)/max(integral_179);
end