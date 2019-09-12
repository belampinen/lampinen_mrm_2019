function imsize = my_imsize(n_cols_figure, height_width_ratio, journal)

if (nargin < 3), journal = []; end
if (nargin < 2), height_width_ratio = []; end

if (isempty(journal)), journal = 'mrm'; end
if (isempty(height_width_ratio)), height_width_ratio = 1; end

% Determine width
switch (lower(journal))
    
    case 'mrm'
        
        switch (n_cols_figure)
            
            case 1
                image_width = 3.42;
                
            case 1.5
                image_width = 5.12;
                
            case 2
                image_width = 6.9;
                
            otherwise
                error('Strange # cols');
        end        
        
    otherwise
        
        error(['Journal not implemented: ' journal]);
end

% Return Size
imsize = [image_width image_width * height_width_ratio];

end