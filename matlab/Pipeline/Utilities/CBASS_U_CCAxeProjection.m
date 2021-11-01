function [db2CoorProj] = CBASS_U_CCAxeProjection(db2Data, db2Centers, blNorm)
% Utilities returning the coordinates of a set of observation contained in
% DB2DATA when projected on the axis linking the centroids contained in
% db2Centers. The coordinates are normalized so that the centroid of H0 and
% Test occupy positon 0 and 1 respectively.

%Deals with optional arguments
narginchk(1, 3)
if nargin < 3, blNorm = false; end

% Checks that the input is properly formatted
[inNCntr, inNPar_Cntr] = size(db2Centers);
if inNCntr < 2, error('There should be at least two centers in db2Center'); end

[inNObs , inNPar] = size(db2Data);
if inNPar_Cntr ~= inNPar
    error('db2Data and db2Centers must have the same number of columns');
end

% Initializes the coordinates matrix
db2CoorProj = nan(inNObs, inNCntr * (inNCntr - 1) / 2);

% Calculates the centroid to centroid coordinates
iCmp = 0;
for iCt1 = 1:inNCntr - 1
    for iCt2 = iCt1 + 1:inNCntr
        iCmp = iCmp + 1;
        db2CoorProj(:, iCmp) = (db2Data - db2Centers(iCt2, :)) * (db2Centers(iCt1, :) - db2Centers(iCt2, :))';
        if blNorm
            db2CoorProj(:, iCmp) = 2 * db2CoorProj(:, iCmp) ./ (norm(db2Centers(iCt1, :) - db2Centers(iCt2, :)).^2) - 1;
        end
    end
end

