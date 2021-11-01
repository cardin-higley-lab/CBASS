function dbPVal = CBASS_U_TestRate(dbRate_1, inN_1, dbRate_2, inN_2)
% Utility testing a difference in binomial rates.

% Checkes input
narginchk(3, 4);
if any([dbRate_1 dbRate_2] < 0) | any([dbRate_1 dbRate_2] > 1)
    error('Rates must be comprised between zero and 1')
end

% If the number of occurence for the sedond rate is not provided, it is
% taken that dbRate_2 is H.0
if nargin < 4
    dbEpsilon = abs(dbRate_1 - dbRate_2)./...
        sqrt(dbRate_2 * (1 - dbRate_2) ./ inN_1);
else % Otherwise compares the two ratios as if they are issued from two differnt samples
    dbRate_Tot = (dbRate_1 * inN_1 + dbRate_2 * inN_2) ./ (inN_1 + inN_2);
    dbEpsilon = abs(dbRate_1 - dbRate_2)./...
        sqrt((dbRate_Tot * (1 - dbRate_Tot) ./ inN_1) + (dbRate_Tot * (1 - dbRate_Tot) ./ inN_2));
end

% Significance of the increase
dbPVal =  2 * (1 - normcdf(dbEpsilon));