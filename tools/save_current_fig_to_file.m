function save_current_fig_to_file (out_name, out_dir, imsize, res, out_format)
% function fig_qute (out_name, out_path, imsize, res)
% By: Filip Szczepankiewicz, 2014.09.24
% This function saves the current graphic in a customizable format!
% out_name      The name of the file that is saved. Do not include file
%               ending!
% out_dir       The path to the directory where the image is saved. If this
%               is undefined, set to '', or left empty the default is PWD.
% imsize        The absolute size [x y] of the image in inches.
% res           The resolution of the image in dots per inch (DPI)
% out_format    The format of the saved image. See code for available
%               formats.


do_fix_axes = 0;
do_disable_col_invert = 1;


if (nargin < 2 || strcmp(out_dir, '') || isempty(out_dir))
    out_dir = '/Volumes/DIFFDATA/BJORN/FEXI/133_FEXI_Isopulse';
end

if (nargin < 3)
    % Keep current size of the figure window
else
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [1 1 imsize]);
end

if (nargin < 4)
    res = 600;
end

if nargin < 5
    out_format = 'tiff';
    
    % Other compatible formats:
    % png, pdf, jpeg, epsc, bmp
end


% Fix missing file separator
if ~strcmp(out_dir(end), '\')
    out_dir = [out_dir '\'];
end

% Fix missing axis
if do_fix_axes
    f = get(gcf);
    for i = 1:numel(f.Children)
        tmp = get(f.Children(i));
        try
            ct = tmp.CameraTarget * (1 + 0.0000001);
            set(f.Children(i), 'CameraTarget', ct)
        catch
            disp('Camera could not be adjusted!')
        end
    end
end

% Enables changed background colors when printing
if do_disable_col_invert
    set(gcf, 'InvertHardCopy', 'off');
end

% Print to file
tmp_nam = [out_dir out_name];

switch version
    case {'8.2.0.701 (R2013b)', '8.0.0.783 (R2012b)'}
        f_num = gcf;
        
    case {'8.4.0.150421 (R2014b)', '8.6.0.267246 (R2015b)'}
        fig = gcf;
        f_num = fig.Number;
end

print(['-d' out_format], ['-f' num2str(f_num)], ['-r' num2str(res)], tmp_nam);



