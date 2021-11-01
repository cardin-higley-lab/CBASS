function CBASS_Plot_PeakHistogram(sPULSE, blLegend)
%Plotting utility of the bout pipeline to plot histogram of the score value
%at peak for the real and surrogate (H0) data issued by CBASS_L3_GetPulse.
%The function will also plot the threshold.

%Handle optional arguments
narginchk(1, 2)
if nargin < 2, blLegend = true; end

%Compute the value of the score at peak for real and surrogate (H0) data
db1PkVal     = sPULSE.db1Score(sPULSE.bl1Peak);
db1PkVal_H0  = sPULSE.db1Score_H0(sPULSE.bl1Peak_H0);

% Calculates the bins of the histogram
blNorm = all([db1PkVal db1PkVal_H0] >= -.5 & [db1PkVal db1PkVal_H0] <= 1);  
if blNorm, db1Bin = -.5:.01:1; else; [~, db1Bin] = hist([db1PkVal db1PkVal_H0], 150); end

%Computes the histograms
db1HistPk = hist(db1PkVal, db1Bin);
db1HistPk_H0 = hist(db1PkVal_H0, db1Bin);

% Plots the histograms
bar(db1Bin, [db1HistPk; db1HistPk_H0]'); hold on
db1YL = ylim; plot([1 1] * sPULSE.dbThrPeak, db1YL, '--r');
if blLegend, legend({'LFP', 'RndLFP', '.05% Threshold'}, 'Location', 'NorthWest'); end
if sPULSE.dbP_KS > 0.0001, chPKS = sprintf('%.4f', sPULSE.dbP_KS); else, chPKS = '<0.0001'; end
dbSigRate = sum(sPULSE.bl1Pulse)./sum(sPULSE.bl1Peak);
if sPULSE.dbP_Rate > 0.0001, chPRate = sprintf('%.4f', sPULSE.dbP_Rate); else, chPRate = '<0.0001'; end
title(sprintf('Peaks: LFP vs Surrogate (KS Test: p = %s;) Pulse Rate: %.2f (p = %s)', chPKS, dbSigRate, chPRate));