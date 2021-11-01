function [cBRANCH, in1NObs_Bch, db1Rate_Bch, db1PVal_Bch] = ...
    GetSigBranch_v2(db2Tree, in1NObs, db1Rate, dbSigThrs)
% Utility meant to be used by CBASS_L2_GetEnrichedRegion. Groups  branches
% if they are significantly enriched as a whole.

% Creates a reference for the indices of the nodes
inNLeaves   = size(db2Tree, 1) + 1;

% Checks that the format of in1NObs and db1Rate is correct
if ~isvector(in1NObs) || numel(in1NObs) ~= inNLeaves, error('in1NObs is improperly formatted'); end 
if ~isvector(db1Rate) || numel(db1Rate) ~= inNLeaves, error('db1Rate is improperly formatted'); end

% Calculates the overall ratio
[~, dbRate_All] = WeightedRate(in1NObs, db1Rate);

% Calculates the significance of each cluster
db1PVal = BinomialTest(in1NObs, db1Rate, dbRate_All);

% Select a starting pool of significant leaves
in1StartLeaves = find(db1PVal < dbSigThrs & db1Rate > dbRate_All);

% Initializes a cell array grouping the clusters belinging to specific
% branches
cBRANCH = {};
[in1NObs_Bch, db1Rate_Bch, db1PVal_Bch] = deal([]);

% Finds branches of the dendrogram where all the leaves are significant
while ~isempty(in1StartLeaves)
    % Initializes the starting leaves and the branch
    in1Leaves       = in1StartLeaves(1);
    inCurrentBranch = in1StartLeaves(1);
    blBchSig        = true;
    
    % Aggregate the branch from the leaves until the whole is not
    % significant
    while blBchSig
        % Finds the connected branch and all its leaves
        [inNode, inLinked]  = FindLinked(db2Tree, inCurrentBranch);
        in1NewLeaves        = FindLeaves([], db2Tree, inLinked);
        
        % Checke if the group is significant as a whole
        [inNObs_Test, dbRate_Test] = WeightedRate(in1NObs([in1Leaves in1NewLeaves]), ...
            db1Rate([in1Leaves in1NewLeaves]));
        if BinomialTest(inNObs_Test, dbRate_Test, dbRate_All) < dbSigThrs & ...
            dbRate_Test > dbRate_All
            % If yes continues
            in1Leaves = cat(2, in1Leaves, in1NewLeaves);
            inCurrentBranch = inNLeaves + inNode; 
        else
            % If not stops
            blBchSig = false;
        end
    end
    % Checks that there leaves do not overlap with existing branches (there
    % should in theory only be one overlapping branch)
    bl1Overlap = cellfun(@(x) any(ismember(in1Leaves, x)), cBRANCH);
    if any(bl1Overlap)
        in1Overlap = find(bl1Overlap);
        cBRANCH{in1Overlap(1)} = unique([cBRANCH{bl1Overlap} in1Leaves]);
        cBRANCH(in1Overlap(2:end)) = [];
    else
        cBRANCH = cat(2, cBRANCH, {in1Leaves});
    end

    % Removes Branch members from the starter leaves
    in1StartLeaves(ismember(in1StartLeaves, in1Leaves)) = [];
end
    
% Calculates branch metrics
for iBch = 1:length(cBRANCH)
    [inNObs_Bch, dbRate_Bch] = WeightedRate(in1NObs(cBRANCH{iBch}), db1Rate(cBRANCH{iBch}));
    dbPVal_Bch = BinomialTest(inNObs_Bch, dbRate_Bch, dbRate_All);
    in1NObs_Bch = cat(2, in1NObs_Bch, inNObs_Bch);
    db1Rate_Bch = cat(2, db1Rate_Bch, dbRate_Bch);
    db1PVal_Bch = cat(2, db1PVal_Bch, dbPVal_Bch);
end

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