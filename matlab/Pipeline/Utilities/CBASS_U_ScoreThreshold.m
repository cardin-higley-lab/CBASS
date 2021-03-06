function dbThreshold = CBASS_U_ScoreThreshold2(db2Data, db1Score)
% CBASS utility. Calculates a threshold for the enrichment score DB1SCORE
% returned by the function CBASS_U_EnrichmentScore. The function returns 
% the threshold that maximizes the normalized mahalanobis distances between 
% the event above and under the threshold. The normalized distance can be 
% thought of as multidimensional analog of the t-statistics. The threshold 
% chosen optimizes the distance while taking sampling variability into 
% account. This version of the function uses fminsearch.
%
% Input -------------------------------------------------------------------
%
% db2Data:      a matrix (observation x parameter) representing the a set
%               of observations in a parameter space
% db1Score:     a score valued between 0 and 1 assigned to each observation

%Performs a search
hFUN        = @(x) - MDNORM(x, db2Data, db1Score);
dbThreshold = fminsearch(hFUN, median(db1Score));

% - Utility function ------------------------------------------------------

function db_MDNorm = MDNORM(dbThr, db2Data, db1Score)
%Sets the distance to infinity if it is out of bound
if dbThr <= 0 || dbThr >= 1, db_MDNorm = -Inf;
else
    %Gets troughs above threshold and their rate
    bl1TSelIdx  = db1Score > dbThr;
    
    %Computes the coordinates of the centroids of trough above and under
    %enrichment score threshold
    db1CntrSup  = mean(db2Data(bl1TSelIdx));
    db1CntrMin  = mean(db2Data(~bl1TSelIdx));
    
    %Computes the euclidean and mahalanobis distance between the centroids
    db_MahalD   = MahalanobisD(db2Data, db1CntrSup - db1CntrMin);
    db_MDNorm   = db_MahalD ./ sqrt(1./sum(bl1TSelIdx) + 1./sum(~bl1TSelIdx));
end

function db1MahaD = MahalanobisD(db2RefPop, db2Obs)
%DB1MAHAD = MahalanobisD(DB2REFPOP, DB2OBS)
%Returns DB1MAHAD: the Mahalanobis distances of a set of observation
%DB2OBS relative to a set of reference observation DB2REFPOP.

db2CCov     = nancov(db2RefPop); %Covariance matrix of the cluster
db1CMu      = nanmean(db2RefPop, 1); %Centroid of the cluster

db2ObsCnt   = db2Obs - db1CMu; %Centered coordintates
db1MahaD    = sqrt(nansum((db2ObsCnt/db2CCov).*db2ObsCnt, 2)); %Squared Mahalonobis distances
