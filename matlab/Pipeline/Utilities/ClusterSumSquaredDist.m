function [db1SSD, db2SSD_Rand] = ClusterSumSquaredDist(db2Data, inNCluMax, inNIter)
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
    
    for iNCl = 2:inNCluMax
        
        % Computes the Log Likelihood of the randomized model
        db2SSD_Rand(iItr, iNCl - 1) = SS_Dist(db2DataRnd, iNCl);
    end
end

% Computes the observation probabilities
for iNCl = 2:inNCluMax
    db1SSD(iNCl - 1) = SS_Dist(db2Data, iNCl);
end

% Calculates the Bayesian information criterion of the clustering
% db1NrmSSDist = (db1SS_D - mean(db2SS_D_Rand))./std(db2SS_D_Rand);


function dbSS_D = SS_Dist(db2Data, iNCl)

% Does the K-Mean correction
[~,  ~, db1SumD]  = kmeans(db2Data, iNCl);
dbSS_D = sum(db1SumD);
