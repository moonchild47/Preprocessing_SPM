function Standard_Preprocessing_Vers2
% SPM12 preprocessing pipeline
% based on scripts by Glad Mihai, Martin Domin,
% Ulrike Horn and Valérie Chanoine
% Date: October 20, 2017
% ulrike.horn@uni-greifswald.de

% general directory structure has to be:
% rootPath
%   -subject01
%       -funcDir
%           -run_001 (session 1)
%           -run_002 (session 2)
%       -structDir
%       -fieldmapDir

% or if there is a fieldmap for every session:
% rootPath
%   -subject01
%       -run_001 (session 1)
%           -funcDir
%           -fieldmapDir
%       -run_002 (session 2)
%           -funcDir
%           -fieldmapDir
%       -structDir

clear
% Initialise SPM
spm('defaults','fmri');
spm_jobman('initcfg');

% Define a general configuration structure 'cfg'
cfg.rootPath        =   'E:\mmInt_Saulin\Data\MRI\SPM';
cfg.funcDir         =   'func';                 % Dir where functional epis are | parent folder: subject (or session)
cfg.batchDir        =   'batch';                % Dir where batch files are saved | will be created if needed
cfg.structDir       =   'mpr_t1_ns_sag_p2_iso';                 % Dir of structural T1 image | parent folder: subject
cfg.fieldMapDir     =   'gre_field_34';             % fieldmap directory | parent folder: subject (or session)
cfg.prefix          =   'func';                    % prefix of all raw data files
cfg.ext             =   '\.*nii';               % extension of files: \.*nii or \.*img
cfg.anatomical      =   'anatomy';              % prefix for the anatomical scan ('s', 'anatomy', ...)

cfg.subjects = {
%     'mm_Int_04'
%     'mm_Int_05'
%     'mm_Int_08'
%     'mm_Int_09'
%     'mm_Int_11'
%     'mm_Int_12' %too many field maps --> took only the first ones (run2/3)
%     'mm_Int_13'
%     'mm_Int_14'
%     'mm_Int_15'
%     'mm_Int_16' %functional directories not correctly named (order ER-MC)
%     'mm_Int_17' %functional directories not correctly named (order ER-MC)
%     'mm_Int_18' %functional directories not correctly named (order ER-MC)
%     'mm_Int_19' %functional directories not correctly named (order ER-MC)
%     'mm_Int_20'
%     'mm_Int_21'
%     'mm_Int_22'
%     'mm_Int_23'
%     'mm_Int_24'
%     'mm_Int_25'
%     'mm_Int_27'
%     'mm_Int_29'
%     'mm_Int_30'
%     'mm_Int_31'
%     'mm_Int_32'
%     'mm_Int_33'
%     'mm_Int_34'
%     'mm_Int_35'
%     'mm_Int_36'
%     'mm_Int_37'
%     'mm_Int_38'
%     'mm_Int_39'
%     'mm_Int_40'
    'mm_Int_41'
    'mm_Int_42'
    'mm_Int_43'
    'mm_Int_44'
    'mm_Int_45'
    'mm_Int_46'
    'mm_Int_47'
    };
cfg.sessions = {
    'bold_MB_Induktion_ER'
    'bold_MB_Verteilung_ER'
    'bold_MB_Induktion_MM'
    'bold_MB_Verteilung_MM'
    };
%% Needed parameters
cfg.multiplefieldmaps   =   0;                      % multiple fieldmaps (one for each session) =1 or only one for all =0
cfg.phaseenc            =   'y';                    % phase encoding direction x, y, z, xy, etc.
cfg.tr                  =   1;                      % repetition time in seconds
cfg.nSlices             =   72;                     % number of slices
cfg.sliceOrder          =   'interleaved';          % descending, ascending, interleaved
cfg.refSlice            =   round(cfg.nSlices/2);   % reference slice for slice timing
cfg.voxStruc            =   [1 1 1];                % voxel size anatomical scans
cfg.voxFunc             =   [2.5 2.5 2.5];                % voxel size functional scans after normalization
cfg.fwhm                =   6;                      % smoothing kernel


%% what do you want to do
% standard procedure: fieldmap, realign, coregister, segment, normalize_own
fieldmap            =   1;  % optional but recommended
realign             =   1;
slicetime           =   0;  % optional
slicetime_before    =   0;  % 1 = slice timing before realignment, 0 = afterwards
coregister          =   1;
segment             =   1;  % new segment (standard in SPM12)
normalize_own       =   0;  % use own template for normalization (only for certain groups as children etc)
normalize_std       =   1;  % normalizes with standard-template from SPM (MNI space)
smooth              =   1;  % do not use when normalizing into your own template space

% If you do not want to use the T1 for normalization use the script
% old_norm_wo_segment.m (not recommended, only for anatomical
% scans with artifacts)

%%
% Loop over subjects
for iSubject=1:numel(cfg.subjects)
    cfg.subPath          = fullfile (cfg.rootPath, cfg.subjects{iSubject});
    cfg.structPath       = fullfile (cfg.subPath, cfg.structDir);
    cfg.batchPath        = fullfile (cfg.subPath, cfg.batchDir);
    if cfg.multiplefieldmaps==1
        funcPath    = fullfile(cfg.subPath, cfg.sessions{1}, cfg.funcDir);
        fieldmapPath = fullfile(cfg.subPath,cfg.sessions{1},cfg.fieldMapDir);
    else
        funcPath    = fullfile(cfg.subPath, cfg.funcDir, cfg.sessions{1});
        fieldmapPath = fullfile(cfg.subPath,cfg.fieldMapDir);
    end
    % check if everything is correct
    if ~exist(funcPath,'dir')
        error('Error. Functional directory not found. Please read comments on top of the script')
    end
    if ~exist(cfg.structPath,'dir')
        error('Error. Anatomical directory not found.')
    end
    if ~exist(fieldmapPath,'dir')&& fieldmap
        error('Error. Fieldmap directory not found.')
    end
    % if there is no batch directory, make one
    type = exist(cfg.batchPath,'dir');
    if type ~= 7
        mkdir(cfg.batchPath);
    end
    %-------------------------------
    % Do Preprocessing step by step
    %-------------------------------
    
    if fieldmap
        DoFieldMap(cfg);
    end
    
    if realign
        if slicetime
            if slicetime_before
                DoSliceTiming(cfg);
                DoRealignUnwarp(cfg);
            else
                DoRealignUnwarp(cfg);
                DoSliceTiming(cfg);
            end
        else
            DoRealignUnwarp(cfg);
        end
    else
        if slicetime
            DoSliceTiming(cfg);
        end
    end
    if segment
        % before doing anything with the anatomical image
        % correct the position so that coregistration afterwards will be
        % fine
        vbm8_set_center_of_mass_nogui(cfg.structPath);
        DoSegment(cfg);
    end
    
    if coregister
        DoCoregister(cfg);
    end
    if normalize_std
        DoNormalize(cfg);
    end
    if smooth
        DoSmooth(cfg);
    end
end
if normalize_own
    DoCreateTemplate(cfg);
    DoDartelNormalize(cfg);
end
end

%%
function DoFieldMap (cfg)
% Parameters used to calculate vdm.
te1         = 4.92;         % short TE
te2         = 7.38;         % long TE
epifm       = 0;          % no epi-based field map, since GRE
tert        = 11.56;       % total epi read out time in ms,
kdir        = -1;          % blip direction (+1/-1), Siemens standard -1
mask        = 0;           % don't mask brain, done by default
match       = 0;          % don't match field map to EPI
write_unw   = 0;      % don't write unwarped epi

%if there are multiple field map files loop over sessions
if cfg.multiplefieldmaps==1
    for isess=1:numel(cfg.sessions)
        fieldmapPath = fullfile(cfg.subPath,cfg.sessions{isess},cfg.fieldMapDir);
        epiPath = fullfile(cfg.subPath,cfg.sessions{isess},cfg.funcDir);
        % Check to see if fieldmap dir exists
        % Skip the vdm calculation if not
        try
            if exist(fieldmapPath,'dir') == 7
                disp('Field Map directory exists.')
                % Check if field map is already calculated. If so, then skip it.
                dirContent = dir(fullfile(fieldmapPath,'vdm*'));
                if isempty(dirContent)
                    fprintf('\nCalculating fieldmap.\nIt takes just a bit...\n')
                    FieldMap_preprocess(fieldmapPath,epiPath,...
                        [te1, te2, epifm, tert, kdir, mask, match, write_unw]);
                else
                    fprintf('\nField Map already calculated; skipping calculation...\n')
                end
            end
        catch err
            disp('=============================================================')
            disp('Cannot find field map directory.')
            disp('  It is either named differently than expected or it doesn''t exist.')
            disp('  The analysis will continue, but check the field map directory before proceeding!')
            disp('=============================================================')
            rethrow(err);
        end
    end
else
    fieldmapPath = fullfile(cfg.subPath,cfg.fieldMapDir);
    epiPath = fullfile(cfg.subPath,cfg.funcDir,cfg.sessions{1});
    % Check to see if fieldmap dir exists
    % Skip the vdm calculation if not
    try
        if exist(fieldmapPath,'dir') == 7
            disp('Field Map directory exists.')
            % Check if field map is already calculated. If so, then skip it.
            dirContent = dir(fullfile(fieldmapPath,'vdm*'));
            if isempty(dirContent)
                fprintf('\nCalculating fieldmap.\nIt takes just a bit...\n')
                FieldMap_preprocess(fieldmapPath,epiPath,...
                    [te1, te2, epifm, tert, kdir, mask, match, write_unw]);
            else
                fprintf('\nField Map already calculated; skipping calculation...\n')
            end
        end
    catch err
        disp('=============================================================')
        disp('Cannot find field map directory.')
        disp('  It is either named differently than expected or it doesn''t exist.')
        disp('  The analysis will continue, but check the field map directory before proceeding!')
        disp('=============================================================')
        rethrow(err);
    end
end
end
%%
function DoRealignUnwarp(cfg)
% assign appropriate vector for the phase encoding direction
if strcmp(cfg.phaseenc,'x')
    wrapping_vec = [1 0 0];
elseif strcmp(cfg.phaseenc,'y')
    wrapping_vec = [0 1 0];
elseif strcmp(cfg.phaseenc,'z')
    wrapping_vec = [0 0 1];
elseif strcmp(cfg.phaseenc,'xy')
    wrapping_vec = [1 1 0];
elseif strcmp(cfg.phaseenc,'yx')
    wrapping_vec = [1 1 0];
elseif strcmp(cfg.phaseenc,'xz')
    wrapping_vec = [1 0 1];
elseif strcmp(cfg.phaseenc,'zx')
    wrapping_vec = [1 0 1];
elseif strcmp(cfg.phaseenc,'yz')
    wrapping_vec = [0 1 1];
elseif strcmp(cfg.phaseenc,'zy')
    wrapping_vec = [0 1 1];
else
    wrapping_vec = [0 0 0];
end

clear matlabbatch
% check if slice timing was done before (are there files with prefix af)
if cfg.multiplefieldmaps==1
    test=fileselector(fullfile(cfg.subPath,...
        cfg.sessions{1},cfg.funcDir),['a',cfg.prefix],cfg.ext);
else
    test=fileselector(fullfile(cfg.subPath,...
        cfg.funcDir,cfg.sessions{1}),['a',cfg.prefix],cfg.ext);
end
if isempty(test)
    prefix = cfg.prefix;
else
    prefix = ['a',cfg.prefix];
end
% loop over sessions to collect all EPIs
for isess=1:numel(cfg.sessions)
    % check if there are multiple fieldmaps
    if cfg.multiplefieldmaps==1
        fieldmapPath = fullfile(cfg.subPath,cfg.sessions{isess},cfg.fieldMapDir);
        funcimages = cellstr(fileselector(fullfile(cfg.subPath,...
            cfg.sessions{isess},cfg.funcDir),prefix,cfg.ext));
    else
        fieldmapPath = fullfile(cfg.subPath,cfg.fieldMapDir);
        funcimages = cellstr(fileselector(fullfile(cfg.subPath,...
            cfg.funcDir,cfg.sessions{isess}),prefix,cfg.ext));
    end
    matlabbatch{1}.spm.spatial.realignunwarp.data(isess).scans = funcimages;
    vdm = fileselector(fieldmapPath,'vdm',cfg.ext);
    % check if the voxel distortion map exists to use it for correction
    if isempty(vdm)
    else
        matlabbatch{1}.spm.spatial.realignunwarp.data(isess).pmscan = cellstr(vdm);
    end
end

matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm = 0;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp = 2;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap = wrapping_vec;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.jm = 0;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.sot = [];
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.rem = 1;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.noi = 5;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap = wrapping_vec;
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.mask = 1;
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';

save(fullfile(cfg.batchPath, 'SPM12_matlabbatch_Realign&Unwarp.mat'),'matlabbatch');
fprintf('\nRunning realign and unwarp\n');
spm_jobman('run',matlabbatch);
end

%%
function DoSliceTiming(cfg)
clear matlabbatch;
% check if realignment has been done before (prefix uf)
% else use raw files (e.g. prefix f)
if cfg.multiplefieldmaps==1
    test = fileselector(fullfile(cfg.subPath,cfg.sessions{1},...
        cfg.funcDir),['u',cfg.prefix],cfg.ext);
else
    test = fileselector(fullfile(cfg.subPath,cfg.funcDir,...
        cfg.sessions{1}),['u',cfg.prefix],cfg.ext);
end
if isempty(test)
    prefix = cfg.prefix;
else
    prefix = ['u',cfg.prefix];
end
% Loop for sessions
matlabbatch{1}.spm.temporal.st.scans = {};
for isess=1:numel(cfg.sessions)
    if cfg.multiplefieldmaps==1
        epi = fileselector(fullfile(cfg.subPath, cfg.sessions{isess},...
            cfg.funcDir),prefix,cfg.ext);
    else
        epi = fileselector(fullfile(cfg.subPath, cfg.funcDir,...
            cfg.sessions{isess}),prefix,cfg.ext);
    end
    matlabbatch{1}.spm.temporal.st.scans{isess} = cellstr(epi);
end
% Slice order
if strcmp(cfg.sliceOrder,'descending')
    matlabbatch{1}.spm.temporal.st.so = cfg.nSlices:-1:1;
elseif strcmp(cfg.sliceOrder,'ascending')
    matlabbatch{1}.spm.temporal.st.so = 1:1:cfg.nSlices;
elseif strcmp(cfg.sliceOrder,'interleaved') && mod(cfg.nSlices,2)==0   % even
    matlabbatch{1}.spm.temporal.st.so = [2:2:cfg.nSlices 1:2:cfg.nSlices];
elseif strcmp(cfg.sliceOrder,'interleaved') && mod(cfg.nSlices,2)~=0   % odd
    matlabbatch{1}.spm.temporal.st.so = [1:2:cfg.nSlices 2:2:cfg.nSlices];
end

matlabbatch{1}.spm.temporal.st.refslice = cfg.refSlice;
matlabbatch{1}.spm.temporal.st.nslices = cfg.nSlices;
matlabbatch{1}.spm.temporal.st.tr = cfg.tr;
matlabbatch{1}.spm.temporal.st.ta = cfg.tr-(cfg.tr/cfg.nSlices);
matlabbatch{1}.spm.temporal.st.refslice = cfg.refSlice;
matlabbatch{1}.spm.temporal.st.prefix = 'a';

save(fullfile(cfg.batchPath, 'SPM12_matlabbatch_SliceTiming.mat'),'matlabbatch');
fprintf('\nRunning Slice Timing\n');
spm_jobman('run',matlabbatch);
end
%%
function DoSegment(cfg)
clear matlabbatch;
% Select the T1 3D image coregistered
anat = cellstr(fileselector(cfg.structPath,cfg.anatomical,cfg.ext));
% Get template of each cerebral tissue
tmpGM           = {fullfile(spm('Dir'),'tpm', 'TPM.nii,1')};
tmpWM           = {fullfile(spm('Dir'),'tpm', 'TPM.nii,2')};
tmpCSF          = {fullfile(spm('Dir'),'tpm', 'TPM.nii,3')};
tmpBone         = {fullfile(spm('Dir'),'tpm', 'TPM.nii,4')};
tmpSoftTissue   = {fullfile(spm('Dir'),'tpm', 'TPM.nii,5')};
tmpAirBck       = {fullfile(spm('Dir'),'tpm', 'TPM.nii,6')};
matlabbatch{1}.spm.spatial.preproc.channel.vols = anat;
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001; % Bias regularisation light
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;  % 60 mm cutoff
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1]; % save bias field and corrected version
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = tmpGM;
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 1];  %native+Dartel
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = tmpWM;
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 1];  %native+Dartel
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = tmpCSF;
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 1];  %native+Dartel
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = tmpBone;
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];  %native
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = tmpSoftTissue;
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];  %native
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = tmpAirBck;
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];  %none
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];   % Inverse + Forward

save(fullfile(cfg.batchPath, 'SPM12_matlabbatch_Segment.mat'),'matlabbatch');
fprintf('\nRunning segmentation\n');
spm_jobman('run',matlabbatch);

clear matlabbatch;
% Image Calculator - Create brain image (skull-stripped bias corrected)
c1 = cellstr(fileselector(cfg.structPath,'c1',cfg.ext));
c2 = cellstr(fileselector(cfg.structPath,'c2',cfg.ext));
c3 = cellstr(fileselector(cfg.structPath,'c3',cfg.ext));
bias= cellstr(fileselector(cfg.structPath,['m',cfg.anatomical],cfg.ext));

matlabbatch{1}.spm.util.imcalc.input = [c1;c2;c3;bias];
matlabbatch{1}.spm.util.imcalc.output = 'Brain';
matlabbatch{1}.spm.util.imcalc.outdir = {cfg.structPath};
matlabbatch{1}.spm.util.imcalc.expression = '(i1 + i2 + i3) .* i4';
save(fullfile(cfg.batchPath, 'SPM12_matlabbatch_Skull.mat'),'matlabbatch');
fprintf('\nCalculating skull stripped bias corrected brain image \n');
spm_jobman('run',matlabbatch);

end
%%
function DoCoregister(cfg)
% Get skull stripped bias corrected brain image
brain = fileselector(cfg.structPath,'Brain',cfg.ext);

% Get mean realigned and unwarped EPI as well as all other unwarped EPIs
if cfg.multiplefieldmaps==1
    meanUnwarp = fileselector(fullfile(cfg.subPath, cfg.sessions{1}, cfg.funcDir),'mean',cfg.ext);
    if size(meanUnwarp,1)>1
        error('Error. \nThere is more than one mean realigned EPI.')
    else
        test = fileselector(fullfile(cfg.subPath, cfg.sessions{1},cfg.funcDir),...
            ['ua',cfg.prefix],cfg.ext);
    end
else
    meanUnwarp = fileselector(fullfile(cfg.subPath, cfg.funcDir, cfg.sessions{1}),'mean',cfg.ext);
    if size(meanUnwarp,1)>1
        error('Error. \nThere is more than one mean realigned EPI.')
    else
        test = fileselector(fullfile(cfg.subPath, cfg.funcDir,cfg.sessions{1}),...
            ['ua',cfg.prefix],cfg.ext);
    end
end
if isempty(test)
    if cfg.multiplefieldmaps==1
        test=fileselector(fullfile(cfg.subPath, ...
            cfg.sessions{1},cfg.funcDir),['au',cfg.prefix],cfg.ext);
    else
        test=fileselector(fullfile(cfg.subPath, ....
            cfg.funcDir,cfg.sessions{1}),['au',cfg.prefix],cfg.ext);
    end
    if isempty(test)
        prefix = ['u',cfg.prefix];
    else
        prefix = ['au',cfg.prefix];
    end
else
    prefix = ['ua',cfg.prefix];
end

epi = {};
% Loop on sessions to collect sliced EPIs
for isess=1:numel(cfg.sessions)
    if cfg.multiplefieldmaps==1
        f = fileselector(fullfile(cfg.subPath,...
            cfg.sessions{isess},cfg.funcDir),prefix,cfg.ext);
    else
        f = fileselector(fullfile(cfg.subPath,...
            cfg.funcDir,cfg.sessions{isess}),prefix,cfg.ext);
    end
    epi = vertcat(epi, cellstr(f));
end
clear matlabbatch
matlabbatch{1}.spm.spatial.coreg.estimate.ref = cellstr(brain);
matlabbatch{1}.spm.spatial.coreg.estimate.source = cellstr(meanUnwarp);
matlabbatch{1}.spm.spatial.coreg.estimate.other = cellstr(epi);
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 ...
    0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

save(fullfile(cfg.batchPath, 'SPM12_matlabbatch_Coregister.mat'),'matlabbatch');
fprintf('\nRunning coregistration\n')
spm_jobman('run',matlabbatch);

end
%%
function DoNormalize(cfg)

% Get Field Deformation image
forwardDeformation = fileselector(cfg.structPath,['y_',cfg.anatomical],cfg.ext);

% Get skull stripped bias corrected brain image
brain = fileselector(cfg.structPath,'Brain',cfg.ext);

% prefix can be uf, uaf or auf
if cfg.multiplefieldmaps==1
    test = fileselector(fullfile(cfg.subPath, cfg.sessions{1},cfg.funcDir),...
        ['ua',cfg.prefix],cfg.ext);
else
    test = fileselector(fullfile(cfg.subPath, cfg.funcDir,cfg.sessions{1}),...
        ['ua',cfg.prefix],cfg.ext);
end
if isempty(test)
    if cfg.multiplefieldmaps==1
        test=fileselector(fullfile(cfg.subPath, ...
            cfg.sessions{1},cfg.funcDir),['au',cfg.prefix],cfg.ext);
    else
        test=fileselector(fullfile(cfg.subPath, ....
            cfg.funcDir,cfg.sessions{1}),['au',cfg.prefix],cfg.ext);
    end
    if isempty(test)
        prefix = ['u',cfg.prefix];
    else
        prefix = ['au',cfg.prefix];
    end
else
    prefix = ['ua',cfg.prefix];
end

epi = {};
% Loop on sessions to collect sliced EPIs
for isess=1:numel(cfg.sessions)
    if cfg.multiplefieldmaps==1
        f = fileselector(fullfile(cfg.subPath,...
            cfg.sessions{isess},cfg.funcDir),prefix,cfg.ext);
    else
        f = fileselector(fullfile(cfg.subPath,...
            cfg.funcDir,cfg.sessions{isess}),prefix,cfg.ext);
    end
    epi = vertcat(epi, cellstr(f));
end
clear matlabbatch;

matlabbatch{1}.spm.spatial.normalise.write.subj.def = {forwardDeformation};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {brain};
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-90 -126 -72;...
    90 90 108];
%matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = NaN(2,3);
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = cfg.voxStruc;
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
% spm_get_defaults('normalise.write.bb')?
matlabbatch{2}.spm.spatial.normalise.write.subj.def = {forwardDeformation};
matlabbatch{2}.spm.spatial.normalise.write.subj.resample = cellstr(epi);
matlabbatch{2}.spm.spatial.normalise.write.woptions.bb = [-90 -126 -72;...
    90 90 108];
%matlabbatch{2}.spm.spatial.normalise.write.woptions.bb = NaN(2,3);
matlabbatch{2}.spm.spatial.normalise.write.woptions.vox = cfg.voxFunc;
matlabbatch{2}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{2}.spm.spatial.normalise.write.woptions.prefix = 'w';

% save(fullfile(cfg.batchPath, 'SPM12_matlabbatch_Normalize.mat'),'matlabbatch');
save(fullfile(cfg.batchPath, 'SPM12_matlabbatch_Normalize.mat'),'matlabbatch');
spm_jobman('run',matlabbatch);
end
%%
function DoCreateTemplate(cfg)
%% DARTEL
rc1 = cell(numel(cfg.subjects),1);
rc2 = rc1;
% batch is in directory of subject 1
batchPath=fullfile(cfg.rootPath, cfg.subjects{1},cfg.batchDir);
% Go through each subject
counter = 1;
for iSubject = 1:numel(cfg.subjects)
    structPath = fullfile (cfg.rootPath,cfg.subjects{iSubject},cfg.structDir);
    if exist(structPath,'dir') == 7
        % Take the rc1*.nii image and add it to a matrix.
        rc1{counter,1} = fileselector(structPath,'rc1',cfg.ext);
        % Take the rc2*.nii image and add it to a matrix.
        rc2{counter,1} = fileselector(structPath,'rc2',cfg.ext);
        counter = counter + 1;
    else
        fprintf('\nDirectory for %s missing.\n',cfg.subjects{iSubject})
    end
end

clear matlabbatch
matlabbatch{1}.spm.tools.dartel.warp.images = {rc1,rc2};
matlabbatch{1}.spm.tools.dartel.warp.settings.template = 'Template';
matlabbatch{1}.spm.tools.dartel.warp.settings.rform = 0;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).rparam = [4 2 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).K = 0;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).slam = 16;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).rparam = [2 1 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).K = 0;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).slam = 8;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).rparam = [1 0.5 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).K = 1;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).slam = 4;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).rparam = [0.5 0.25 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).K = 2;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).slam = 2;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).K = 4;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).slam = 1;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).K = 6;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).slam = 0.5;
matlabbatch{1}.spm.tools.dartel.warp.settings.optim.lmreg = 0.01;
matlabbatch{1}.spm.tools.dartel.warp.settings.optim.cyc = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.optim.its = 3;

type = exist(batchPath,'dir');
if(type ~= 7)
    mkdir(batchPath);
end
save(fullfile(batchPath,sprintf('SPM12_Create_template_%s.mat',date)),'matlabbatch');
fprintf('\nCreating template.\n')
spm_jobman('run',matlabbatch);
end
%%
function DoDartelNormalize(cfg)
template = fileselector(fullfile(cfg.rootPath,cfg.subjects{1},...
    cfg.structDir),'Template_6',cfg.ext);
for iSubject = 1:numel(cfg.subjects)
    clear matlabbatch
    subPath          = fullfile (cfg.rootPath, cfg.subjects{iSubject});
    structPath       = fullfile (subPath, cfg.structDir);
    batchPath        = fullfile (subPath, cfg.batchDir);
    if exist(subPath,'dir') == 7
        % prefix can be uf, uaf or auf
        if cfg.multiplefieldmaps==1
            test = fileselector(fullfile(subPath, ...
                cfg.sessions{1},cfg.funcDir),['ua',cfg.prefix],cfg.ext);
        else
            test = fileselector(fullfile(subPath, ...
                cfg.funcDir,cfg.sessions{1}),['ua',cfg.prefix],cfg.ext);
        end
        if isempty(test)
            if cfg.multiplefieldmaps==1
                test=fileselector(fullfile(subPath, ...
                    cfg.sessions{1},cfg.funcDir),['au',cfg.prefix],cfg.ext);
            else
                test=fileselector(fullfile(subPath, ...
                    cfg.funcDir,cfg.sessions{1}),['au',cfg.prefix],cfg.ext);
            end
            if isempty(test)
                prefix = ['u',cfg.prefix];
            else
                prefix = ['au',cfg.prefix];
            end
        else
            prefix = ['ua',cfg.prefix];
        end
        epi = {};
        % Loop on sessions to collect unwarped EPIs
        for isess=1:numel(cfg.sessions)
            if cfg.multiplefieldmaps==1
                f = fileselector(fullfile(subPath,...
                    cfg.sessions{isess},cfg.funcDir),prefix,cfg.ext);
            else
                f = fileselector(fullfile(subPath,...
                    cfg.funcDir,cfg.sessions{isess}),prefix,cfg.ext);
            end
            epi = vertcat(epi, cellstr(f));
        end
        flowfield = fileselector(structPath,'u_rc1',cfg.ext);
        
        % normalize all functional images
        matlabbatch{1}.spm.tools.dartel.mni_norm.template = {template};
        matlabbatch{1}.spm.tools.dartel.mni_norm.data.subj.flowfield = {flowfield};
        matlabbatch{1}.spm.tools.dartel.mni_norm.data.subj.images = cellstr(epi);
        matlabbatch{1}.spm.tools.dartel.mni_norm.vox = cfg.voxFunc;
        matlabbatch{1}.spm.tools.dartel.mni_norm.bb = [NaN NaN NaN
            NaN NaN NaN];
        % spm_get_defaults('normalise.write.bb')?
        matlabbatch{1}.spm.tools.dartel.mni_norm.preserve = 0;
        matlabbatch{1}.spm.tools.dartel.mni_norm.fwhm = ones(1,3)*cfg.fwhm;
        
        % normalize skull stripped bias corrected brain image
        brain = fileselector(cfg.structPath,'Brain',cfg.ext);
        %         anat_image = fileselector(structPath,['m',cfg.anatomical],cfg.ext);
        matlabbatch{2}.spm.tools.dartel.mni_norm.template = {template};
        matlabbatch{2}.spm.tools.dartel.mni_norm.data.subj.flowfield = {flowfield};
        matlabbatch{2}.spm.tools.dartel.mni_norm.data.subj.images = cellstr(brain);
        matlabbatch{2}.spm.tools.dartel.mni_norm.vox = [NaN NaN NaN];
        matlabbatch{2}.spm.tools.dartel.mni_norm.bb = [NaN NaN NaN
            NaN NaN NaN];
        matlabbatch{2}.spm.tools.dartel.mni_norm.preserve = 0;
        matlabbatch{2}.spm.tools.dartel.mni_norm.fwhm = [1 1 1];
        
        % save batch
        type = exist(batchPath,'dir');
        if type ~= 7
            mkdir(batchPath);
        end
        save(fullfile(batchPath,sprintf('SPM12_Norm2MNI_%s_%s.mat',date)),'matlabbatch');
        fprintf('\nNormalising subject %s to MNI',cfg.subjects{iSubject})
        spm_jobman('run',matlabbatch);
        
        % rename the functional images to see which smoothing kernel was
        % applied
        fprintf('\nRenaming files\n')
        for isess=1:numel(cfg.sessions)
            if cfg.multiplefieldmaps==1
                epiPath = fullfile(subPath,cfg.sessions{isess},cfg.funcDir);
            else
                epiPath = fullfile(subPath,cfg.funcDir,cfg.sessions{isess});
            end
            clear matlabbatch
            matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {epiPath};
            matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '^sw';
            matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';
            matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_move.files(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^sw)',...
                substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
            matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.moveto = {epiPath};
            matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.pattern = 'sw';
            matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.patrep.repl = sprintf('s%02dw',cfg.fwhm);
            matlabbatch{2}.cfg_basicio.file_dir.file_ops.file_move.action.moveren.unique = false;
            
            spm_jobman('run',matlabbatch);
        end
    else
        sprintf('\nDirectory %s for %s missing.\n',cfg.subjects{iSubject})
    end
end

end
%%
function DoSmooth(cfg)
clear matlabbatch;
P = {};
% Get Normalized EPI files of all sessions
% for isess=1:numel(cfg.sessions)
%     if cfg.multiplefieldmaps==1
%         f = fileselector(fullfile(cfg.subPath,...
%             cfg.sessions{isess},cfg.funcDir),'w',cfg.ext);
%     else
%         f = fileselector(fullfile(cfg.subPath,...
%             cfg.funcDir,cfg.sessions{isess}),'w',cfg.ext);
%     end
%     P = vertcat(P,cellstr(f));
% end
for isess=1:numel(cfg.sessions)
    if cfg.multiplefieldmaps==1
        f = fileselector(fullfile(cfg.subPath,...
            cfg.sessions{isess},cfg.funcDir),'w',cfg.ext);
    else
        f = fileselector(fullfile(cfg.subPath,...
            cfg.funcDir,cfg.sessions{isess}),'w',cfg.ext);
    end
    P = vertcat(P,cellstr(f));
end
matlabbatch{1}.spm.spatial.smooth.data = cellstr(P);
matlabbatch{1}.spm.spatial.smooth.fwhm = ones(1,3)*cfg.fwhm;
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = sprintf('s%02d',cfg.fwhm);

% save(fullfile(cfg.batchPath, 'SPM12_matlabbatch_Smooth.mat'),'matlabbatch');
save(fullfile(cfg.batchPath, 'SPM12_matlabbatch_Smooth.mat'),'matlabbatch');
spm_jobman('run',matlabbatch);
end