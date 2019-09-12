function smr_optimize_log(fid, s_message, do_disp, do_add_newline, do_nothing, do_add_timestamp)
% function my_log(fid, s_message, do_disp, do_add_newline, do_nothing, do_add_timestamp)
%
% Logs a message to (an already open) file and (potentially) displays it 
% on screen


if (nargin < 2), fid = -1; end
if (nargin < 3), do_disp = 0; end
if (nargin < 4), do_add_newline = 1; end
if (nargin < 5), do_nothing = 0; end
if (nargin < 6), do_add_timestamp = 1; end

% do_nothing is n/verbose_level, so that if it is above one, do nothing)
if (do_nothing > 1), return; end


if (do_add_newline)
    s_newline = char([13 10]);
else
    s_newline = [];
end

if (do_add_timestamp)
    s_timestamp = [datestr(now, 'YYYY-mm-DD HH:MM') ' '];
else
    s_timestamp = [];
end

if (fid ~= -1)
    fwrite(fid, [s_timestamp s_message s_newline]);
end

if (do_disp), disp(s_message); end
end