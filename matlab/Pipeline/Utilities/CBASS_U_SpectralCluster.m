function [in1Member, db2Vec, db1Val] = CBASS_U_SpectralCluster(db2Adj, chOpt)
%Utility performing spectral clustering on the similarity matrix DB2SIM

%Calculates the size of db2Sim
[inN, inNCol] = size(db2Adj);
if inN ~= inNCol || ~ismatrix(db2Adj), error('db2Sim must be a square matrix'); end

%Checks for optional argument
narginchk(1, 2);
if nargin < 2; chOpt = 'unorm'; end
if ~ismember(chOpt, {'unorm', 'norm', 'rw'})
    fprintf('option string non recognized treated as unormalized');
end

%Calculates the degree matrix
db2Deg  = diag(sum(db2Adj, 2));

%Calculates the Laplacian
db2Lap = db2Deg - db2Adj;
switch chOpt
    case 'norm' % normalized Laplacian
        db2Lap = sqrt(db2Deg) \ db2Lap / sqrt(db2Deg);
    case 'rw' % random walk Laplacian
        db2Lap = db2Deg \ db2Lap;
end

% Performs the eigenvalue decomposition of the Laplacian
[db2Vec,db2Val] = eig(db2Lap, eye(length(db2Lap)), 'chol');
db2Val = real(db2Val);
% Sorts the values of the Laplacian
[db1Val, in1Srt]  = sort(diag(db2Val)');
db2Vec = db2Vec(:, in1Srt);

% Gets the number of clusters
in1Zero     = find(db1Val < 10^-10);
in1Member   = kmeans(db2Vec(:, in1Zero), length(in1Zero));