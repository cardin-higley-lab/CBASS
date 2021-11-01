function in1Score = CBASS_U_EnrichmentScore(db2Data, bl1State, bl1Baseline, inNClu, dbSigThrs, blZScore, inNIter)
% CBASS utility. Calculates an enrichment score for each observations (i.e.
% rows) in the data matrix DB2DATA. The score represent how likely an
% observation is to fall into a region that has more observation labelled
% by the boolean BL1STATE than by chance. The score is estimated by
% performing INNITR random partition of the data, and testing for a higher
% rate of occurence of BL1STATE within each region using a binomial test.
% Partitions are performed using a method analogous to the first step of
% the k-means algorithm. INNCLU centers are drawn from the observation at
% random and observations are assigned to their closest center.
%
% Input -------------------------------------------------------------------
%
% db2Data:      a matrix (observation x parameter) representing the a set
%               of observations in a parameter space
% bl1State:     a boolean vector indexing each observation
% bl1Baseline:  (optional) a boolean vector indexing each observation and
%               having no overlap with bl1State. Default is ~bl1State.
% inNClu:       (optional) the number of regions used to conmpute the
%               score. Default is 20. Higher values will yield steeper
%               score distributions.
% dbSigThrs:    (optional) threshold for the significance of the enrichment
%               in each region.  P-Values are computed with a binomial 
%               test. Default is 10.^-4.
% bl1Zscore:    (optional) logical specifying if data is to be zscored for 
%               partitioning. Default is true.
% inNItr:       (optional) integer specifying how many iteration of
%               the algorithm are performed. Default is 1000;

% Checks input
narginchk(2, 7)
if ~exist('bl1Baseline', 'var'), bl1Baseline = ~bl1Epoch; elseif isempty(bl1Baseline), bl1Baseline = ~bl1Epoch; end
if ~exist('blZScore', 'var'), blZScore = true; elseif isempty(blZScore), blZScore = true; end
if ~exist('inNClu', 'var'), inNClu = 20; elseif isempty(inNClu), inNClu = 20; end
if ~exist('dbSigThrs', 'var'), dbSigThrs = 10.^-4; elseif isempty(dbSigThrs), dbSigThrs = 10.^-4; end
if ~exist('inNMaxIter', 'var'), inNIter = 1000; elseif isempty(inNIter), inNIter = 1000; end

% Checks argurments
inNObs = size(db2Data, 1);
if ~isvector(bl1State) || numel(bl1State) ~= inNObs
    error('bl1State must be a boolean vector having as many elements as there is rows in db2Data')
end
if isrow(bl1State), bl1State = bl1State'; end
if ~isvector(bl1Baseline) || numel(bl1Baseline) ~= inNObs
    warning('bl1Baseline must be a boolean vector the same size aas bl1State. Set to default');
    bl1Baseline = ~bl1State;
end
if isrow(bl1Baseline), bl1Baseline = bl1Baseline'; end
if any(bl1Baseline & bl1State)
    warning('Bl1Baseline and bl1State cannot have common elements. Set to default');
    bl1Baseline = ~bl1State;
end

% Calculates the global rate
dbRate_All      = sum(bl1State)/sum(bl1State | bl1Baseline);

% Zscores data if needed
if blZScore, db2Data = zscore(db2Data); end

% Initializes the count vector
in1SigCount = zeros(inNObs, 1);

% For a number of time set by inNItr:
% 1. Draw inNClu centers from the distribution
% 2. Assign each points to its closest center to define a region
% 3. Checks whether each region is enriched
% 4. Adds a count for each point situated in an enriched region
for iItr = 1:inNIter
    % Select intialization centroids
    db2Cntr = DrawCentroids(db2Data, inNClu);
    
    % Calculates the distance of each point to the centroid and assigns each
    % point to its cluster
    db2D           = pdist2(db2Data, db2Cntr);
    [~, in1Clu]    = min(db2D, [], 2);
    
    % Calculates the enrichment of each regions and its significance
    for iClu = 1:inNClu
        bl1Clu = in1Clu == iClu;
        inNObs = sum(bl1Clu);
        dbRate = sum(bl1State(bl1Clu))/sum(bl1State(bl1Clu) | bl1Baseline(bl1Clu));
        if dbRate > dbRate_All & BinomialTest(inNObs, dbRate, dbRate_All) < dbSigThrs
            in1SigCount(bl1Clu) = in1SigCount(bl1Clu) + 1;
        end
    end
end

% Normalizes the count by the number of interation to obtain an enrichment
% score for each point
in1Score = in1SigCount ./ inNIter;

% Utilities ---------------------------------------------------------------

% Function 1 - Pick a set of centroids spannin the volume of the space
function db2Cntr = DrawCentroids(db2Data, inNClu)

% Select inNCLU random observation to serve as centers
in1Idx  = datasample(1:size(db2Data, 1), inNClu, 'Replace', false);
db2Cntr = db2Data(in1Idx, :);

% Function 2 - Performs a Binomial test for one sample and a test rate
function dbXPVal = BinomialTest(inXNObs, dbXRate, dbXRate_Tot)
dbXPVal = 2 .* (1 - normcdf(abs(dbXRate - dbXRate_Tot)./...
    sqrt(dbXRate_Tot .* (1 - dbXRate_Tot) ./ inXNObs)));