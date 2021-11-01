function [status, commandOut] = CBASS_L2_PlotEmbed(sTROUGH, chTmpPath, chOutPath, chExperiment, chLabelTag, ...
    chMethod, chDataFormat, blZScore, inN_Component, db1Label, chFormatImg, chRotate3D, chAddLegend, inFontSize, chDiscrete)
% Wrapper function for matlab calling the python function CBASS_L2_PlotEmbed.py
% We use this function for to plot the 
%
% Input -------------------------------------------------------------------
%
% sTROUGH:          a structure containing the following fields:
%                   -.db1FilterBand an (1 x 2) array describing the frequency
%                   band of interest i.e. [30 80] for 30 to 80 Hz.
%                   -.db2Trough  a (2 * channel x trough) matrix containing the
%                   hilbert transform of each channel of sREC.db2LFP filtered
%                   in the band defined in sTROUGH.db1FilterBand at the trough
%                   of the filtered signal in the reference channel inRefChan
%                   -.in1Index the indices of the trough in sREC.db2LFP
% chTmpPath:        a character array indicating the temporary directory where
%                   files used by the embedding function will be saved.
% chOutPath:        a character array indicating the path to the output folder, 
%                   in which the embedded data will be saved.
% chExperiment:     character array indicating the name of the experiment and any
%                   other information that should be contained in the name of 
%                   saved file.
% chLabelTag        character array dedicated to labelling embedding plots using
%                   different plotting settings.
% chMethod:         character array specifying the method used for embedding. 
%                   The options are UMAP, PHATE and PCA. The default is 'UMAP'.      
% chDataFormat:     (optional) a character array specifying the format of
%                   the hilbert transforms output. Can be 'complex' or
%                   'polar'. Default is 'complex'.
% blZScore:         logical specifying if trought data is to be zscored for 
%                   k-means partitioning. Default is True.
% inN_Component:    (optional) a number between one and three specifying the 
%                   dimensionality of the embedded space. The default is a '3'.
% db1Labels:        array containing the labels associated to the data samples being 
%                   plotted. It is used to defined the color code of the plot. 
%                   The default is coloring samples according to their odering 
%                   (eg, timestamps)
% chFormatImg:      character array specifying the format of the output image.
%                   The options are 'png', 'eps, and 'gif'. The default is 'png'.
% chRotate3D:       character array specifying if the 3D plot will be saved as rotating 
%                   gif image or single image (png or eps). The default is 'False'.
% chAddLegend:      character array to hide ('False') or include ('True') the legend
%                   to the figure. The default is 'True'.
% inFontSize:       integer specifying fontsize of the plot. The default is 10.
% chDiscrete:       string specifying the settings for the legend. Use 'True' for discrete colorbar 
%                   and 'False' is for continuous. The default is 'True'.

% Output ------------------------------------------------------------------
% status            command exit status. When the command is successful, status is 0. 
%                   Otherwise, status is a nonzero integer.
% commandOut        Output of the operating system command, returned as a character vector.

% Checks argument
narginchk(5, 15)
if ~exist('chMethod', 'var'),       chMethod = 'umap';      elseif isempty(chMethod),       chMethod = 'umap'; end
if ~exist('chDataFormat', 'var'),   chDataFormat = 'polar'; elseif isempty(chDataFormat),   chDataFormat = 'polar'; end
if ~exist('blZScore', 'var'),       blZScore = true;        elseif isempty(blZScore),       blZScore = true; end
if ~exist('inN_Component', 'var'),  inN_Component = 3;      elseif isempty(inN_Component),  inN_Component = 3; end
if ~exist('chFormatImg', 'var'),    chFormatImg = 'png';    elseif isempty(chFormatImg),    chFormatImg = 'png'; end
if ~exist('chRotate3D', 'var'),     chRotate3D = 'False';   elseif isempty(chRotate3D),     chRotate3D = 'False'; end
if ~exist('chAddLegend', 'var'),    chAddLegend = 'True';   elseif isempty(chAddLegend),    chAddLegend = 'True'; end
if ~exist('inFontSize', 'var'),     inFontSize = 10;        elseif isempty(inFontSize),     inFontSize = 10; end
if ~exist('chDiscrete', 'var'),     chDiscrete = 'True';    elseif isempty(chDiscrete),     chDiscrete = 'True'; end
if ~exist('db1Label', 'var'),       db1Label = []; end

% Checks the input
if ~ismember(chMethod, {'umap', 'phate', 'pca'})
    fprintf('%s is not a valid method set to umap', chMethod); chMethod = 'umap'; 
end
if ~ismember(chDataFormat, {'polar', 'complex'})
    fprintf('%s is not a valid method set to umap', chDataFormat); 
    chDataFormat = 'complex'; 
end

% Export variables to temp for test
data            = sTROUGH.db2Trough;
    % Convert data to the wanted representation if wanted (the complex
    % trough can be represented in complex or polar coordinates)
if strcmp(sTROUGH.chDataFormat, 'complex') && strcmp(chDataFormat, 'polar')
    inNChan = size(data, 2)/2;
    data    = [abs(data(:, 1:inNChan) + 1i * data(:, inNChan + 1: end)), ...
        angle(data(:, 1:inNChan) + 1i * data(:, inNChan + 1: end))]; 
elseif strcmp(sTROUGH.chDataFormat, 'polar') && strcmp(chDataFormat, 'complex')
    inNChan     = size(data, 2)/2;
    db2DatCmplx = data(:, 1:inNChan) .* exp(1i * data(:, inNChan + 1:end));
    data        = [real(db2DatCmplx) imag(db2DatCmplx)]; 
end
if blZScore, data = zscore(data); end %zscore the data if required
time_stamps     = sTROUGH.in1Index - 1; %Time stamps of interest adjusted for python indexing
% Sets the labels
if isempty(db1Label), labels = 1:size(data,1); % Just the frames for now
else, labels = db1Label; end
method_embed    = chMethod; % or 'umap' default: phate
n_components    = inN_Component;
% chFormatImg     = 'png'
% chRotate3D        = 'False'
% chAddLegend      = 'False'
category        = -1; %To select running, for example, "-1" will just use frame numbers for labeling

save(fullfile(chTmpPath,['tmp_vars_' chExperiment '_' chLabelTag '.mat']), 'data', 'time_stamps', 'labels')

% Run this if you just want to plot the embedding
tic
disp(['Plotting the embedded data...'])
commandStr = "python ../Pipeline/CBASS_L2_PlotEmbed.py " + ... 
    " --tmpPath " + chTmpPath + ...
    " --experiment " + chExperiment + ...
    " --n_components " + num2str(n_components) + ...
    " --method " + method_embed + ...
    " --category " + num2str(category) + ...
    " --label_name " + chLabelTag + ...
    " --add_legend " + chAddLegend + ...
    " --format_img " + chFormatImg + ...
    " --fontsize " + inFontSize + ...
    " --discrete " + chDiscrete + ...
    " --rotate3D " + chRotate3D + ...
    " --outPath " + chOutPath + "" ;
% [status, commandOut] = system(sprintf(commandStr));
[status, commandOut] = system(sprintf('conda activate gammaBouts_env && %s', commandStr));
commandStr
if status ~= 0 % For some reason, the cluster only works with 'source'. I am adding this as an option.
    [status, commandOut] = system(sprintf('source activate gammaBouts_env && %s', commandStr));
end
toc