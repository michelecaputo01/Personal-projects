function x_3D = reshape_into_3D(x, lines, columns, DEBUG)

% Convert a 2D matrix x (with shape depth x lines*columns) to a 3D tensor 
% (with shape lines x columns x depth). If DEBUG == 1, the correctness of
% the conversion is also checked

if (size(x,2) ~= lines*columns)
    error('x, lines, columns variables are incompatible')
end


% Perform reshape

depth = size(x,1);
x_3D  = reshape(x', [lines columns depth]);


% Check reshape

if (DEBUG)
    
    mismatches = 0;
    index = 1;
    
    for (j = 1:columns)
        for (i = 1:lines)
            for (k = 1:depth)
                
                if(x(k,index) ~= x_3D(i,j,k))
                    mismatches = mismatches+1;
                end
                
            end
            
            index = index+1;
        end
    end
    
    if (mismatches)
        error(['Passing from image_2D to image_3D, there are ', num2str(mismatches), ' mismatches'])
    end
end