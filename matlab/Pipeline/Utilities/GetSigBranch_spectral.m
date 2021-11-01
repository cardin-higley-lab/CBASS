function [cBRANCH, in1NObs_Bch, db1Rate_Bch, db1PVal_Bch] = ...
    GetSigBranch_spectral(db2Data, in1NObs, db1Rate, dbSigThrs)

% Creates a reference for the indices of the nodes
inNLeaves   = size(db2Data, 1);

% Checks that the format of in1NObs and db1Rate is correct
if ~isvector(in1NObs) || numel(in1NObs) ~= inNLeaves, error('in1NObs is improperly formatted'); end 
if ~isvector(db1Rate) || numel(db1Rate) ~= inNLeaves, error('db1Rate is improperly formatted'); end

% Calculates the overall ratio
[~, dbRate_All] = WeightedRate(in1NObs, db1Rate);

% Calculates the significance of each cluster
db1PVal = BinomialTest(in1NObs, db1Rate, dbRate_All);

% Select a starting pool of significant leaves
in1StartLeaves = find(db1PVal < dbSigThrs & db1Rate > dbRate_All);

%First assume all the cluster are significant, but check for how many very
%small eigenvalues we have. That should be the real number of clusters
% [in1Clu,V,D] = spectralcluster(db2Data(in1StartLeaves,:),length(in1StartLeaves), 'Distance','precomputed','LaplacianNormalization','symmetric');
[in1Clu,V,D] = spectralcluster(db2Data(:,2:end),size(db2Data,1), 'Distance','precomputed','LaplacianNormalization','symmetric');
real_num_of_clusters = find(D < 10e-5);

% [in1Clu,V,D] = spectralcluster(db2Data(in1StartLeaves,2:end),length(real_num_of_clusters), 'Distance','precomputed','LaplacianNormalization','symmetric');
[in1Clu,V,D] = spectralcluster(db2Data(:,2:end),length(real_num_of_clusters), 'Distance','precomputed','LaplacianNormalization','symmetric');


% Only keep the enriched regions
in1Clu = in1Clu(in1StartLeaves);
% Initializes a cell array grouping the clusters belinging to specific
% branches
cBRANCH = {};
[in1NObs_Bch, db1Rate_Bch, db1PVal_Bch] = deal([]);

% Creates a cell array of branches
inNClu=max(in1Clu);
cBRANCH = cell(1, inNClu);
for iBch = 1:inNClu, cBRANCH{iBch} = find(in1Clu == iBch); end

% Calculates branch metrics
[in1NObs_Bch, db1Rate_Bch, db1PVal_Bch] = deal(nan(size(cBRANCH)));
for iBch = 1:inNClu
    [in1NObs_Bch(iBch), db1Rate_Bch(iBch)] = WeightedRate(in1NObs(cBRANCH{iBch}), db1Rate(in1StartLeaves(cBRANCH{iBch})));
    db1PVal_Bch(iBch) = BinomialTest(in1NObs_Bch(iBch), db1Rate_Bch(iBch), dbRate_All);
end

% Only keep significant branches
bl1Keep = db1PVal_Bch < dbSigThrs & db1Rate_Bch > dbRate_All;
cBRANCH     = cBRANCH(bl1Keep);
in1NObs_Bch = in1NObs_Bch(bl1Keep);
db1Rate_Bch = db1Rate_Bch(bl1Keep);
db1PVal_Bch = db1PVal_Bch(bl1Keep);

% Sorts the branches by overal hit rate
[db1Rate_Bch, in1SrtIdx] = sort(db1Rate_Bch, 'descend');
cBRANCH     = cBRANCH(in1SrtIdx);
in1NObs_Bch = in1NObs_Bch(in1SrtIdx);
db1PVal_Bch = db1PVal_Bch(in1SrtIdx);

function [inNObs, dbRatio] = WeightedRate(in1NObs, db1Rate)
% Calculates the total number of observation INNOBS and the ratio of hits
% DBRATIO resulting in grouping clusters having numbers of observations
% and hit ratio stored in IN1NOBS and DB1RATIO respectively
inNObs  = sum(in1NObs);
dbRatio = sum(in1NObs .* db1Rate) ./ inNObs;

function db1PVal = BinomialTest(in1NObs, db1Rate, dbRate_Tot)
% Return the p-values of a binomial test  of the difference between
% clusters having hit rate DBRATIO
db1PVal = 2 * (1 - normcdf(abs(db1Rate - dbRate_Tot)./...
            sqrt(dbRate_Tot * (1 - dbRate_Tot) ./ in1NObs)));

function [inNode, inLinked] = FindLinked(db2Tree, inCurrentBranch)
inNode = find(db2Tree(:, 1) == inCurrentBranch);
if ~isempty(inNode), inLinked = db2Tree(inNode, 2); return; end
inNode = find(db2Tree(:, 2) == inCurrentBranch);
inLinked = db2Tree(inNode, 1);

function in1Leaves  = FindLeaves(in1Leaves, db2Tree, inStartNode)
if inStartNode > size(db2Tree, 1) + 1
    in1Leaves  = FindLeaves(in1Leaves, db2Tree, db2Tree(inStartNode - size(db2Tree, 1) - 1, 1));
    in1Leaves  = FindLeaves(in1Leaves, db2Tree, db2Tree(inStartNode - size(db2Tree, 1) - 1, 2));
else
    in1Leaves = cat(2, in1Leaves, inStartNode);
end