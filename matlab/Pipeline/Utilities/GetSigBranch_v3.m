function [cBRANCH, in1NObs_Bch, db1Rate_Bch, db1PVal_Bch] = ...
    GetSigBranch_v3(db2Data, in1NObs, db1Rate, dbSigThrs)
% Utility meant to be used by CBASS_L2_GetEnrichedRegion. Sets a global
% threshold based on the ratio of the within and between distance to
% identify branches.

% Creates a reference for the indices of the nodes
inNLeaves   = size(db2Data, 1);

% Conputes the tredd
db2Tree = linkage(db2Data, 'ward');

% Checks that the format of in1NObs and db1Rate is correct
if ~isvector(in1NObs) || numel(in1NObs) ~= inNLeaves, error('in1NObs is improperly formatted'); end 
if ~isvector(db1Rate) || numel(db1Rate) ~= inNLeaves, error('db1Rate is improperly formatted'); end

% Calculates the overall ratio
[~, dbRate_All] = WeightedRate(in1NObs, db1Rate);

% Gets the number of branches (i.e. clusters)
[db1SSD, db2SSD_Rand] = WardSumSquaredDist(db2Data);
db1Nrm_SSD      = (db1SSD - mean(db2SSD_Rand))./std(db2SSD_Rand);
[~, inMinIdx]   = min(db1Nrm_SSD);
inNClu          = inMinIdx + 1;
in1Clu          = cluster(db2Tree, 'maxclust', inNClu);

% Creates a cell array of branches
cBRANCH = cell(1, inNClu);
for iBch = 1:inNClu, cBRANCH{iBch} = find(in1Clu == iBch); end

% Calculates branch metrics
[in1NObs_Bch, db1Rate_Bch, db1PVal_Bch] = deal(nan(size(cBRANCH)));
for iBch = 1:inNClu
    [in1NObs_Bch(iBch), db1Rate_Bch(iBch)] = WeightedRate(in1NObs(cBRANCH{iBch}), db1Rate(cBRANCH{iBch}));
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