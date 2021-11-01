function [db1SSD, db2SSD_Rand] = WardSumSquaredDist(db2Data, inNCluMax, inNIter)
% Utility returning the normalized log likelihood of a clustering partition.
% The partition is compared to the average log likelihood of a partition on
% scrambled data

% Controls for optional arguments
narginchk(1, 3);
if nargin < 2, inNCluMax = round(size(db2Data, 1)./3); end
if nargin < 3, inNIter = 1000; end

% Gets the size of the data
[inNObs, inNPar] = size(db2Data);

% Initializes the random cost matrix
db1SSD         = nan(1, inNCluMax - 1); 
db2SSD_Rand    = nan(inNIter, inNCluMax - 1);

for iItr = 1:inNIter
    % Randomized the data and computes a dendrogram
    db2DataRnd =  db2Data;
    for iPar = 1:inNPar
        db2DataRnd(:, iPar) = db2DataRnd(randperm(inNObs), iPar);
    end
    
    % Compute the tree for the randomized data
    db2Tree = linkage(db2DataRnd, 'ward');
    
    % Calculates the sum squared distance on the randomized data
    for iNCl = 2:inNCluMax
        % Computes the sum squared distance for the randomized model
        db2SSD_Rand(iItr, iNCl - 1) = SS_Dist(db2DataRnd, db2Tree, iNCl);
    end
end

% Computes the ree for the real data
 db2Tree = linkage(db2DataRnd, 'ward');

% Computes the sum squared distance on the real data
for iNCl = 2:inNCluMax
    db1SSD(iNCl - 1) = SS_Dist(db2Data, db2Tree, iNCl);
end

% Calculates the Bayesian information criterion of the clustering
% db1NrmSSDist = (db1SS_D - mean(db2SS_D_Rand))./std(db2SS_D_Rand);


function dbSS_D = SS_Dist(db2Data, db2Tree, inNCl)

% Does computes the sum squared distance
in1CluWard  = cluster(db2Tree, 'maxclust', inNCl);
db1SD   = nan(size(db2Data, 1), 1);
for iCl = 1:inNCl
    bl1Clu          = in1CluWard == iCl;
    db1Ctr          = mean(db2Data(bl1Clu, :));
    db1SD(bl1Clu)   = sum((db2Data(bl1Clu, :) - db1Ctr) .^ 2, 2); 
end
dbSS_D = sum(db1SD);
