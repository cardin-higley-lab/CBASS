function [cBRANCH, in1NObs_Bch, db1Rate_Bch, db1PVal_Bch] = ...
    GetSigBranch(db2Tree, in1NObs, db1Rate, dbSigThrs)
% Utility meant to be used by CBASS_L2_GetEnrichedRegion. Groups branches if
% all the leaves of that branch are significantly enriched.

% Creates a reference for the indices of the nodes
inNLeaves   = size(db2Tree, 1) + 1;

% Calculates the overall ratio
[~, dbRate_All] = WeightedRate(in1NObs, db1Rate);

% Checks the significance of individual nodes
db1PVal = BinomialTest(in1NObs, db1Rate, dbRate_All); 
in1Sig  = find(db1Rate > dbRate_All & db1PVal < dbSigThrs);

% Initializes a cell array grouping the clusters belinging to specific
% branches
cBRANCH = {};
% Finds branches of the dendrogram where all the leaves are significant
while ~isempty(in1Sig)
    in1Leaves       = in1Sig(1);
    inCurrentBranch = in1Sig(1);
    blAllSig        = true;
    while blAllSig
        [inNode, inLinked]  = FindLinked(db2Tree, inCurrentBranch);
        in1NewLeaves        = FindLeaves([], db2Tree, inLinked);
        if all(ismember(in1NewLeaves, in1Sig))
            in1Leaves = cat(2, in1Leaves, in1NewLeaves);
            inCurrentBranch = inNLeaves + inNode; 
        else
            blAllSig = false;
        end
    end
    cBRANCH = cat(2, cBRANCH, {in1Leaves});
    in1Sig(ismember(in1Sig, cBRANCH{end})) = [];
end

% Calculates branch metrics
[in1NObs_Bch, db1Rate_Bch, db1PVal_Bch] = deal(nan(size(cBRANCH)));
for iBch = 1:length(cBRANCH)
    [in1NObs_Bch(iBch), db1Rate_Bch(iBch)] = WeightedRate(in1NObs(cBRANCH{iBch}), db1Rate(cBRANCH{iBch}));
    db1PVal_Bch(iBch) = BinomialTest(in1NObs_Bch(iBch), db1Rate_Bch(iBch), dbRate_All);
end

% Sorts the branches by overal hit rate
[db1Rate_Bch, in1SrtIdx] = sort(db1Rate_Bch, 'descend');
cBRANCH     = cBRANCH(in1SrtIdx);
in1NObs_Bch = in1NObs_Bch(in1SrtIdx);
db1PVal_Bch = db1PVal_Bch(in1SrtIdx);

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