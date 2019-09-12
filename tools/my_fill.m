function str_out = my_fill(str_in, n, do_shrink, my_char, do_center) % add spaces to achieve lenght n

if (nargin < 5), do_center = []; end
if (nargin < 4), my_char = []; end
if (nargin < 3), do_shrink = []; end

if (isempty(do_center)), do_center = 0; end
if (isempty(my_char)), my_char = ' '; end
if (isempty(do_shrink)), do_shrink = 1; end

    
n_in = numel(str_in);
if (n_in < n || do_shrink)
    
    str_out = char(zeros(1,n)) + char(my_char);
    start_index = 1;
    stop_index  = min(n_in, n);
        
    if (do_center && n_in < n)
        start_index = ceil((n - n_in) / 2);
        stop_index  = start_index + n_in - 1; 
    end
    
    if (n_in > 0)
        str_out(start_index:stop_index) = str_in;
    end
    str_out = char(str_out(1:n));
    
else
    
    str_out = str_in;
    
end
end