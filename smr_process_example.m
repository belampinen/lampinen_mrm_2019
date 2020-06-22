function smr_process_example
%
% Purpose: 
%   Process raw (NIFTI) data into SMR maps. This example is based on the 
%   'in vivo' protocol presented in Lampinen et al (2020) MRM
%
% Assumptions: 
%   1) Data consists of one nii file for each combination of 
%      b_delta and TE (here four)
%   2) A corresponding _gdir.txt file exists in the same folder for each such file
%   3) A reversed-polarity counterpart exists for one file
%
% Steps:    
%   1) Identify files and merge them into a single nii file together with a 
%   single xps file.
%   2) Extrapolation-based motion correction (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0141825)
%   3) TOPUP using the reversed-polarity counterpart
%   4) GIBBS correction
%   5) SMR fitting (includes masking and smoothing)
%

%%% ------------------------------------------------------------------------
% PREPARE 
%%% ------------------------------------------------------------------------


%%% User provided details [CHANGE THIS]
%
data_dir   = 'path_to_data';
output_dir = 'path_to_output';
%
file_signature_list = {...                  % strings uniquely identifying each NII file
    'FWFSMS_TE63.nii', ...      
    'FWFSMS_TE130.nii', ...     
    'FWFSMS_TE85.nii', ...      
    'FWFSMS_TE85_PrTE.nii', ... 
    };
b_delta     = [1.0  1.0     1.0 0.6];       % b_delta value for each file
te          = [63   130     85  85]*1e-3;   % TE value for each file
%
pa_signature = 'FWFSMS_TE63_rev_ph';        % reversed-pl
te_pa       = 63*1e-3;


%%% Default
%
opt                 = mdm_opt;
opt.do_overwrite    = 0;
opt.verbose         = 1;
dmri_tag            = 'dmri';


%%% ------------------------------------------------------------------------
% PERFORM 
%%% ------------------------------------------------------------------------

if (~exist(output_dir, 'dir')), mkdir(output_dir); end


%%% 1) Merge files into signle nii file and its corresponding xps file
%
% Identify filenames
d = dir(data_dir);
fn_list = cell(size(file_signature_list));
for c_sign = 1:numel(file_signature_list)
    for c = 1:numel(d)
        if (~isempty(strfind(d(c).name, file_signature_list{c_sign})))
            fn_list{c_sign} = ...
                fullfile(data_dir, d(c).name);
            break;
        end
    end
    if (isempty(fn_list{c_sign}))
        error(['Could not find: ' file_signature_list{c_sign}]);
    end
end
%
% Merge the nii files and create a merged xps
clear s x;
for c_f = 1:numel(fn_list)
    s.nii_fn = fn_list{c_f};
    s.xps       = mdm_xps_from_gdir(mdm_fn_nii2gdir(s.nii_fn), [], b_delta(c_f));
    s.xps.te    = zeros(size(s.xps.b)) + te(c_f);
    x{c_f}      = s;
end
s = mdm_s_merge(x, output_dir, dmri_tag, opt);
s.xps = rmfield(s.xps, 'bt2');              % needs to be removed for the topup step to work


%%% 2) Extrapolation-based motion correction
%
% Subsample LTE with lower b
s_sub = mdm_s_subsample(s, (s.xps.b_delta > 0.9) & (s.xps.s_ind == 1) & ...
    (s.xps.b < 1.1e9), output_dir, opt);
%
% Reference, using the subsample
p = elastix_p_affine(50);
p_fn = elastix_p_write(p, 'p.txt');
s_ref = mdm_mec_b0(s_sub, p_fn, output_dir, opt);
%
% EBMC for the full data
s = mdm_mec_eb(s, s_ref, p_fn, output_dir, opt);


%%% 3) TOPUP (AP against PA)
%
revpol_fn       = msf_find_fn(data_dir, ['*' pa_signature '.nii.gz']);
revpol_xps      = mdm_xps_from_gdir(mdm_fn_nii2gdir(revpol_fn));
revpol_xps.te   = zeros(size(revpol_xps.b)) + te_pa;
mdm_xps_save(revpol_xps, strrep(revpol_fn, '.nii.gz', '_xps.mat'), opt);
s_revpol        = mdm_s_from_nii(revpol_fn);
topup_fn        = fullfile(output_dir, 'dmri_mc_topup.nii.gz');
s_topup         = mdm_s_topup(s, s_revpol, output_dir, topup_fn, opt);


%%% 4) GIBBS
%
gibbs_fn = fullfile(output_dir, 'dmri_mc_topup_gibbs.nii.gz');
if (exist(gibbs_fn, 'file') && ~opt.do_overwrite)
    disp(['Skipping, output file already exists: ' gibbs_fn]);
else
    [I, h] = mdm_nii_read(topup_fn);
    I_corr = mio_gibbs_correction(I);
    mdm_nii_write(I_corr, gibbs_fn, h);
    mdm_xps_save(s_topup.xps, strrep(gibbs_fn, '.nii.gz', '_xps.mat'), opt);
end


%%% 5) Fit SMR
%
dtd_smr_pipe(s, output_dir, opt);
disp('Finished fitting SMR');


end