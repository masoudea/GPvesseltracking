function[direction] = check_direction(direction)

    if direction > 180
        direction= direction-360;
    elseif direction <= -180
        direction= direction+360;
    end
    
end