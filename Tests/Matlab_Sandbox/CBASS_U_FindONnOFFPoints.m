function [in1OnPts, in1OffPts] = CBASS_U_FindONnOFFPoints(bl1EpochVector)
%Synopsis: [IN1ONPTS, IN1OFFPTS] = CBASS_U_FindONnOFFPoints(BL1EPOCHVECTOR)
%Returns the first (IN1ONPTS) and last (IN1OFFPTS) indices of the epochs
%where the boolean vector BL1EPOCHVECTOR is true;
%
% %2016-09-27 QP: Created

% Checks that the input vector is properly formated
if ~isvector(bl1EpochVector), error('The input is not vector'); end
if any(~ismember(unique(bl1EpochVector), [0 1])), error('The input vector must only consist of ones and zeros'); end

% Makes sure the input is a row vector
bl1EpochVector = bl1EpochVector(:);

% Find On indices
in1OnPts    = 1 + find(diff(bl1EpochVector) == 1);
in1OffPts   = find(diff(bl1EpochVector) == -1);

% Checks the ends of bl1EpochVector are not themselves on and off points
if bl1EpochVector(1) == 1; in1OnPts = [1; in1OnPts]; end
if bl1EpochVector(end) == 1; in1OffPts = [in1OffPts; length(bl1EpochVector)]; end