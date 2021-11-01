function [sREGION, sCLU, sFILTERS, cBRANCH] = CBASS_L2_GetEnrichedRegion(sTROUGH, db2LFP, inSampleRate, bl1Epoch, ...
    blZScore, blUseRate, inMethod, inNClu, dbSigThrs, inNMaxIter, blVerbose)
% LEGACY L2 of the bout pipeline: Generate template motifs of activity
% in the band of interst enriched for the state indexed by the logical
% vector bl1Epoch. The algorithm separates a set of events of the input
% signal sREC.db2LFP into an aribitrary number of clusters using the
% k-means algorithm. These events are the output of CBASS_L1_GetTrough and
% correspond to the trougths of oscilatory activity at the band of interest
% in a reference channel. They are represented as the Hilbert transform of
% the band filtered activity of each channel. Clusters are then sorted as a
% function of the fraction of the events they comprise that occured during
% the state defined by bl1Epoch. Clusters that are significantly 'enriched'
% in events occuring during this state are grouped into regions based on
% their adjacency by one of three dendrogram based methods (set by the
% variable inMethod -- see below). Demdrogram are obtained with Ward's
% method.
%
% NOTE: DEPRECATED FUNCTION - This function is an early version of the L2
% of the pipeline kept for legacy. It will give results that are less
% robust and less interpretable than CBASS_L2_GetFilters.m. We advise
% against using it.
%
% Input -------------------------------------------------------------------
%
% sTROUGH:      the output of CBASS_L1_GetTrough (i.e.) a structure
%               requiring the following fields:
%               -.db1FilterBand an (1 x 2) array describing the frequency
%               band of interest i.e. [30 80] for 30 to 80 Hz.
%               -.db2Trough  a (2 * channel x trough) matrix containing the
%               hilbert transform of each channel of sREC.db2LFP filtered
%               in the band defined in sTROUGH.db1FilterBand at the trough
%               of the filtered signal in a reference channel (see
%               CBASS_L1_GetTrough for more detail)
%               -.in1Index the indices of the trough in sREC.db2LFP
% bl1Epopch:    a logical vector, containing as many elements as time
%               samples in db2LFP, indexing the state in which enriched
%               band specific activity is observed.
% blZScore:     (optional) logical specifying if trought data is to be
%               zscored for k-means partitioning and Ward's clustering.
%               Default is true.
% blUseRate:    (optional) logical specifying whether the rate of
%               of enrichment of k-means clusters should be used for Ward's
%               clustering. Default is false.
% inMethod:     (optional) single integer with value 1, 2 or 3. Determines
%               the method used to group k-means cluster into regions. All
%               methods start with building a dendrogram of the k-means
%               cluster using Ward's method. The input is the centroid of
%               each clusters in the space of the input matrix
%               sTROUGH.db2Trough and, optionally, the ratio of events
%               occuring during the state indexed by bl1Epoch:
%               --> Method 1 groups branches if all the leaves of that
%               branch are significantly enriched.
%               --> Method 2 groups branches if they are significantly
%               enriched as a whole.
%               --> Method 3 sets a global threshold based on the ratio of
%               the within and between distance to identify branches.
%               All 3 methods tend to yield several regions which will in
%               turn tend to yield similar filters and overlapping sets of
%               pulses. Default is 2.
% inNClu:       (optional) number of cluster used for k-means partitioning.
%               Default is 20
% dbSigThrs:    (optional) threshold for the significance of the enrichment
%               in trough partition.  P-Values are computed with a binomial 
%               test. Default is 10.^-4.
% inNMaxIter:   (optional) maximum iteration used by the k-means algorithm.
%               Default is 10000.
% blVerbose:    (optional) logical setting whether processing updates
%               should be displayed in the command window. Default is true.
%
% Output ------------------------------------------------------------------
%
% sREGION     	a structure array storing data regarding the enriched
%               regions of the trough data sets and having the following
%               fields:
%               -.bl1Member a boolean indexing troughs used to build the
%               template motif
%               -.inNObs the number of troughs in sFILTER.bl1Member.
%               -.dbRate the fraction of the troughs occuring during the
%               state of interest.
%               -.dbRate_Dev the enrichment for the state of interest i.e.
%               sFILTER.dbRate minus the overall ratio of troughs occuring
%               during the state of interest
%               -.dbPVal the p-value of a binomial test of the of dbRate
%               compared to the overall rate of occurence of the state of
%               interest.
% sCLU          a structure array storing data regarding the clusters
%               generated by the k-means partitioning and having the
%               following fields:
%               -.bl1Member a boolean indexing troughs used to build the
%               template motif
%               -.inNObs the number of troughs in sCLU.bl1Member.
%               -.dbRate the fraction of the troughs occuring during the
%               state of interest.
%               -.dbRate_Dev the enrichment for the state of interest i.e.
%               sCLU.dbRate minus the overall ratio of troughs occuring
%               during the state of interest
%               -.dbPVal the p-value of a binomial test of the of dbRate
%               compared to the overall rate of occurence of the state of
%               interest.

% Sets optional parameters if not provided
narginchk(2, 9);
if ~exist('blZScore', 'var'), blZScore = true; elseif isempty(blZScore), blZScore = true; end
if ~exist('blUseRate', 'var'), blUseRate = true; elseif isempty(blUseRate), blUseRate = false; end
if ~exist('inMethod', 'var'), inMethod = 2; elseif isempty(inMethod), inMethod = 2; end
if ~exist('inNClu', 'var'), inNClu = 20; elseif isempty(inNClu), inNClu = 20; end
if ~exist('dbSigThrs', 'var'), dbSigThrs = 10.^-4; elseif isempty(dbSigThrs), dbSigThrs = 10.^-4; end
if ~exist('inNMaxIter', 'var'), inNMaxIter = 10000; elseif isempty(inNMaxIter), inNMaxIter = 10000; end
if ~exist('blVerbose', 'var'), blVerbose = true; elseif isempty(blVerbose), blVerbose = true; end

% Makes sure that utilities are added to the paths
addpath('./Utilities')

% initializes the random number generator
rng(1949);  
 
% Formats the data
if blZScore, db2Data = zscore(sTROUGH.db2Trough);
else, db2Data = sTROUGH.db2Trough; end

% Print informations on the screen
if blVerbose
    fprintf('\r------ %d CLUSTERS ---------------------------\r', inNClu)
end

% Cluster the data
[in1CluKM, db2Ctr] = kmeans(db2Data, inNClu, 'MaxIter', inNMaxIter);

%%% Checks if clusters are enriched for running the data
% Computes the indices of the bouts where the mouse was running
bl1T_Epoch  = bl1Epoch(sTROUGH.in1Index);

% Computes the overall ratio of running points in the data
dbRate_All      = mean(bl1T_Epoch);

% Initialize a cluster structure containing stats for enrichment
sCLU = struct('bl1Member', cell(1, inNClu), ...
    'inNObs', cell(1, inNClu), ...
    'dbRate', cell(1, inNClu), ...
    'dbRate_Dev', cell(1, inNClu), ...
    'dbPVal', cell(1, inNClu));

sFILTERS = {};

% Loops over clusters
for iClu = 1:inNClu
    bl1Member   = in1CluKM == iClu; % indices of the cluster
    inNObs      = sum(bl1Member); % number of points in the cluster
    sCLU(iClu).bl1Member    = bl1Member;
    sCLU(iClu).inNObs       = inNObs;
    
    % Calculate the ratio of running points for the cluster and the p-value
    % of that ratio under H0 i.e. the hypothesis that the ratio is the same
    % that in the overall population. The p-value is calculated using a
    % binomial distribution.
    dbRate                  = mean(bl1T_Epoch(logical(bl1Member))); % Calculates the ratio
    dbRate_Dev              = dbRate - dbRate_All;
    sCLU(iClu).dbRate       = dbRate; % Stores the ratio
    sCLU(iClu).dbRate_Dev   = dbRate_Dev; % Stores the ratio
    sCLU(iClu).dbPVal       = nanmin(1, 2 * (1 - normcdf(abs(dbRate_Dev)./...
        sqrt(dbRate_All * (1 - dbRate_All) ./ inNObs))));
end

% Calculates the FDR corrected p-values (%WARNING: requires fdr_bh from the
% matlab file exchange)
cSIG = num2cell(fdr_bh([sCLU.dbPVal]));
[sCLU.blSig_FDR] = cSIG{:};

%%% Print the results
% Sort the cluster by ratio
if blVerbose
    [~, in1Sort] = sort([sCLU.dbRate_Dev]);
    
    % Initialize a significance statement array
    cSTATEMENT = {'NON SIG' 'SIG'};
    
    % Loops throught clusters
    for iClu = in1Sort
        fprintf('Cluster %2d (N = %d):\tRate: %.4f\tDev: %.3f\t(p = %.4f)\t%s\r', ...
            iClu, sCLU(iClu).inNObs, sCLU(iClu).dbRate, sCLU(iClu).dbRate_Dev, sCLU(iClu).dbPVal, ...
            cSTATEMENT{sCLU(iClu).blSig_FDR + 1});
    end
end

%%% Links the centroid by distances and identify enriched branches
% Creates a tree of the data
% if blUseRate, db2MatWard = [db2Ctr, [sCLU.dbRate]', ConfusionMatrix]; else, db2MatWard = db2Ctr; end
if blUseRate, db2MatWard = [[sCLU.dbRate]', ConfusionMatrix]; else, db2MatWard = db2Ctr; end
if blZScore, db2MatWard = zscore(db2MatWard); end

% Gets the branches of the dendrogram that are significant
in1NObs = [sCLU.inNObs]; db1Rate = [sCLU.dbRate];
% inMethod = 4;
switch inMethod
    case 1
        % --> Method 1 groups branches if all the leaves of that branch are 
        % significantly enriched 
        db2Tree     = linkage(db2MatWard, 'ward');
        [cBRANCH, in1NObs_Bch, db1Rate_Bch, db1PVal_Bch] = ...
            GetSigBranch(db2Tree, in1NObs, db1Rate, dbSigThrs);
    case 2
        % --> Method 2 groups branches if they are significantly enriched 
        % as a whole 
        db2Tree     = linkage(db2MatWard, 'ward');
        [cBRANCH, in1NObs_Bch, db1Rate_Bch, db1PVal_Bch] = ...
            GetSigBranch_v2(db2Tree, in1NObs, db1Rate, dbSigThrs);
    case 3
        % --> Method 3 sets a global threshold based on the ratio of the
        % within to between distance to identify branches.
        [cBRANCH, in1NObs_Bch, db1Rate_Bch, db1PVal_Bch] = ...
            GetSigBranch_v3(db2MatWard, in1NObs, db1Rate, dbSigThrs);
    case 4
        [cBRANCH, in1NObs_Bch, db1Rate_Bch, db1PVal_Bch] = ...
            GetSigBranch_spectral(db2MatWard, in1NObs, db1Rate, dbSigThrs);
        
end
% Initialize a region structure containing stats for enrichment
inNBch  = length(cBRANCH);
sREGION = struct('bl1Member', cell(1, inNBch), ...
    'inNObs', cell(1, inNBch), ...
    'dbRate', cell(1, inNBch), ...
    'dbRate_Dev', cell(1, inNBch), ...
    'dbPVal', cell(1, inNBch));

for iBch = 1:inNBch
    bl1Member   = ismember(in1CluKM, cBRANCH{iBch}); % indices of the cluster
    sREGION(iBch).bl1Member     = bl1Member; % Membership boolean vector
    sREGION(iBch).inNObs        = in1NObs_Bch(iBch); 
    sREGION(iBch).dbRate        = db1Rate_Bch(iBch);
    sREGION(iBch).dbRate_Dev    = db1Rate_Bch(iBch) - dbRate_All;
    sREGION(iBch).dbPVal        = db1PVal_Bch(iBch);
end

if blVerbose
    fprintf('\r------ %d REGIONS ---------------------------\r', inNBch)
    
        % Loops throught clusters
    for iBch = 1:inNBch
        fprintf('Region %2d (N = %d):\tRate: %.4f\tDev: %.3f\t(p = %.4f)\r', ...
            iBch, sREGION(iBch).inNObs, sREGION(iBch).dbRate, sREGION(iBch).dbRate_Dev, sREGION(iBch).dbPVal);
    end
end

function [inNObs, dbRatio] = WeightedRate(in1NObs, db1Rate)
% Calculates the total number of observation INNOBS and the ratio of hits
% DBRATIO resulting in grouping clusters having numbers of observations
% and hit ratio stored in IN1NOBS and DB1RATIO respectively
inNObs  = sum(in1NObs);
dbRatio = sum(in1NObs .* db1Rate) ./ inNObs;
