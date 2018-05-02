function [seedpoint_position,initial_direction,initial_diameter] = seed_point_selection(input_image)
%     seed_point_selection plots input image and allows us to choose 
%     initial seed-points manually. The function returns initial 
%     seed-point, direction, and diameter.

    figure(1); imshow(uint8(input_image))
    seedpoint = ginput(2);
    hold on
    figure(1);plot(seedpoint(:,1),seedpoint(:,2),'ro')
    drawnow

    dir = atand((-seedpoint(2,2)+ seedpoint(1,2)) / (seedpoint(2,1) - seedpoint(1,1)));
    if seedpoint(1,1) > seedpoint(2,1)
        if seedpoint(1,2) < seedpoint(2,2)
            dir = dir - 180;
        end
        if seedpoint(1,2) > seedpoint(2,2)
            dir = dir + 180;
        end
        if dir ==0
            if seedpoint(1,1) > seedpoint(2,1)
                dir = 180;end

            if seedpoint(1,1) < seedpoint(2,1)
                dir = 0;end
        end
    end
    initial_direction = dir;
    seedpoint_position = [seedpoint(1,1) seedpoint(1,2)];

    [y_cordin(1) x_cordin(1)] =(ginput(1));
    figure(1);plot(y_cordin,x_cordin,'g.')
    [y_cordin(2) x_cordin(2)] =(ginput(1));
    figure(1);plot(y_cordin,x_cordin,'g.')
    initial_diameter = sqrt((y_cordin(1) - y_cordin(2))^2 + (x_cordin(1) - x_cordin(2))^2);

end