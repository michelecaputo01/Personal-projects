function x_2D = reshape_into_2D(x, DEBUG)

% Perform reshape

[lines, columns, depth] = size(x);

x_2D = reshape(x, [(lines*columns) depth])';

if(DEBUG)
    
    x_DEBUG = reshape_into_3D(x_2D, lines, columns, 0);
    
    if(any(any(any(x ~= x_DEBUG))))
        error('The reshape was not performed correctly');
    end
end