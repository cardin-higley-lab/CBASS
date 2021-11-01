function db1Score = CBASS_U_MultiChannelTemplateMatching_2(db2Signal, db2Template, blCenter, blNormalize)
%Synopsis: DB1SCORE = CBASS_U_MultiChannelTemplateMatching_2(DB2SIGNAL, DB2TEMPLATE, [BLNORM])
%Returns a score DB1SCORE indicative of how well the multi-channel signal
%DB2SIGNAL matches the spatio temporal template DB2TEMPLATE
%Input: -DB2SIGNAL a (channel x time sample) matrix
%       -DB2TEMPLATE a (channel x time sample) matrix. The number of
%       channels must be the same as in DB2SIGNAL. The number of time
%       sample must be inferior to half of the duration of the signal.
%       -BLCENTER substract the minimum of the signal defualt is true
%       -BLNORM normalizes the score by norm of db2Signal - default is true
%Output:-DBSCORE a (1 x time sample) row vector. The score represents how
%       well DB2SIGNAL matches DB2TEMPLATE around each time sample and is
%       the dot product S . T where T represents is linearized template
%       (DB2TEMPLATE(:) and S rempresents the linearized chunk of DB2SIGNAL
%       matching the size of DB2TEMPLATE centered on each time points. T is
%       normalized so that all its elements sum to 1.

% Checks the arguments
narginchk(2, 4)
if nargin < 3; blCenter = true; end
if nargin < 4; blNormalize = true; end

% Checks the format of the input
if ~ismatrix(db2Signal), error('db2Signal  must be two dimensional'); end
if ~ismatrix(db2Template), error('db2Template  must be two dimensional'); end
[inNChan, inNSamp] = size(db2Signal);
if size(db2Template, 1) ~= inNChan, error('db2Signal and db2Template must have the same number of channels'); end
inNSmpTmp = size(db2Template, 2);
if xor(isreal(db2Signal), isreal(db2Template)) 
    error('db2Template and db2Signal have to be both real or both complex'); 
end
if inNSmpTmp > inNSamp/2; error('db2Signal must have at least twice as many samples as db2Template'); end
if numel(blCenter) ~= 1; error('blNormalize must a single element boolean'); end
if numel(blNormalize) ~= 1; error('blNormalize must a single element boolean'); end

% If the inputs are complex rearrange them so that the real and imaginary
% parts are each treated as a channel
if ~isreal(db2Signal)
    db2Signal   = [real(db2Signal); imag(db2Signal)];
    db2Template = [real(db2Template); imag(db2Template)];
    inNChan     = inNChan * 2;
end

% Pads db2Signal with zeros for ease of computation if the number of
% sample
% of the template is uneven padds in an asymetric way
inPadBeg = floor(inNSmpTmp/2); 
inPadEnd = ceil(inNSmpTmp/2);
db2Signal = [zeros(inNChan, inPadBeg) db2Signal zeros(inNChan, inPadEnd)];

% Normalizes and transposes the template
if blCenter, db2Template = db2Template - min(db2Template(:)); end
db2Template = db2Template' ./ norm(db2Template(:));

% Computes the centering offset if needed
if blCenter
    dbXCntr = min(db2Signal(:, 1:end - inNSmpTmp));
    for iSmp = 2:inNSmpTmp
        dbXCntr = min(dbXCntr, min(db2Signal(:, iSmp:end - inNSmpTmp + iSmp - 1)));
    end
else
    dbXCntr = 0;
end

% Computes the score
db1Score = db2Template(1, :) * (db2Signal(:, 1:end - inNSmpTmp) - dbXCntr);
if blNormalize, db1SS = sum((db2Signal(:, 1: end - inNSmpTmp) - dbXCntr) .^ 2); end
for iSmp = 2:inNSmpTmp
   db1Score = db1Score + (db2Template(iSmp, :) * (db2Signal(:, iSmp:end - inNSmpTmp + iSmp - 1) - dbXCntr));
   if blNormalize, db1SS = db1SS + sum((db2Signal(:, iSmp: end - inNSmpTmp + iSmp - 1) - dbXCntr) .^ 2); end
end
if blNormalize, db1Score = db1Score ./ sqrt(db1SS); end