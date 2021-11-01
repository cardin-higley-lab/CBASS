function [status, commandOut] = CBASS_L2_Embed(sTROUGH, chTmpPath, chOutPath, chExperiment, chMethod, chDataFormat, blZScore, inN_Component) 
% Wrapper function for matlab calling the python function CBASS_L2_Embed.py
% To visualize the output of each step of our enriched band specific activity motif detection, 
% we project the data to a lower dimension space via dimensionality % reduction. This process 
% is dependent on the integration with Python described in our tutorial 
% (https://github.com/ahof1704/gamma_bouts/wiki/Tutorial#3-python-integration-for-embedding-and-plots)
% 
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
% chMethod:         character array specifying the method used for embedding. 
%                   The options are UMAP, PHATE and PCA. The default is 'UMAP'.      
% blZScore:         logical specifying if trought data is to be zscored for 
%                   k-means partitioning. Default is True.
% chDataFormat:     (optional) a character array specifying the format of
%                   the hilbert transforms output. Can be 'complex' or
%                   'polar'. Default is 'complex'.
% inN_Component:    (optional) a number between one and three specifying the 
%                   dimensionality of the embedded space. The default is a '3'.
%
% Output ------------------------------------------------------------------
% status            command exit status. When the command is successful, status is 0. 
%                   Otherwise, status is a nonzero integer.
% commandOut        Output of the operating system command, returned as a character vector.


% Checks argument
narginchk(4, 8)
if ~exist('chMethod', 'var'),       chMethod = 'umap';          elseif isempty(chMethod),       chMethod = 'umap'; end
if ~exist('chDataFormat', 'var'),   chDataFormat = 'complex';   elseif isempty(chDataFormat),   chDataFormat = 'complex'; end
if ~exist('blZScore', 'var'),       blZScore = true;            elseif isempty(blZScore),       blZScore = true; end
if ~exist('inN_Component', 'var'),  inN_Component = 3;          elseif isempty(inN_Component),  inN_Component = 3; end

% Checks the input
if ~ismember(chMethod, {'umap', 'phate','pca'})
    fprintf('%s is not a valid method set to umap', chMethod); chMethod = 'umap'; 
end
if ~ismember(chDataFormat, {'polar', 'complex'})
    fprintf('%s is not a valid method set to umap', chDataFormat); 
    chDataFormat = 'complex'; 
end

% Export variables to temp for test
data            = sTROUGH(1).db2Trough;
    % Convert data to the wanted representation if wanted (the complex
    % trough can be represented in complex or polar coordinates)
if strcmp(sTROUGH(1).chDataFormat, 'complex') && strcmp(chDataFormat, 'polar')
    inNChan = size(data, 2)/2;
    data    = [abs(data(:, 1:inNChan) + 1i * data(:, inNChan + 1: end)), ...
        angle(data(:, 1:inNChan) + 1i * data(:, inNChan + 1: end))]; 
elseif strcmp(sTROUGH(1).chDataFormat, 'polar') && strcmp(chDataFormat, 'complex')
    inNChan     = size(data, 2)/2;
    db2DatCmplx = data(:, 1:inNChan) .* exp(1i * data(:, inNChan + 1:end));
    data        = [real(db2DatCmplx) imag(db2DatCmplx)]; 
end
if blZScore, data = zscore(data); end %zscore the data if required
time_stamps     = sTROUGH(1).in1Index - 1; %Time stamps of interest adjusted for python indexing
labels          = 1:size(data,1); % Just the frames for now
method_embed    = chMethod; % or 'umap' default: phate
n_components    = inN_Component;
category        = -1; %To select running, for example, "-1" will just use frame numbers for labeling

save(fullfile(chTmpPath,['tmp_vars_' chExperiment '.mat']), 'data', 'time_stamps','labels')

% Get embedding
tic
disp(['Embedding data with ' method_embed ' ...'])
commandStr = "python ../Pipeline/CBASS_L2_Embed.py " + ... 
    " --tmpPath " + chTmpPath + ...
    " --experiment " + chExperiment + ...
    " --n_components " + num2str(n_components) + ...
    " --method " + method_embed + ...
    " --outPath " + chOutPath + "" ;
commandStr
[status, commandOut] = system(sprintf('conda activate gammaBouts_env && %s', commandStr));
if status ~= 0 % For some reason, the cluster only works with 'source'. I am adding this as an option.
    [status, commandOut] = system(sprintf('source activate gammaBouts_env && %s', commandStr));
end
toc