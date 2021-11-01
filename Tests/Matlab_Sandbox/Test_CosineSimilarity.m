% Test -- 1. Cosine similiarity vs normalization
% Set up the cosine similarity function
hCOS = @(x) x' * x ./ (sqrt(diag(x'*x)) .* sqrt(diag(x'*x))');

% Calculates the cosine similarities between a set of random vectors
db2Vec  = rand(100, 1000) + (rand(1, 1000) .* 100 - 50);
cCOS{1} = hCOS(db2Vec);
cCOS_LBL{1} = 'Raw';

% Normalizes the vector by its max and min and calcuates the cosine
% similiarities
db2Vec_MnMx = (db2Vec - min(db2Vec)) ./ (max(db2Vec) - min(db2Vec));
cCOS{2}     = hCOS(db2Vec_MnMx);
cCOS_LBL{2} = 'MinMax';

% Centers the vector by its mean
db2Vec_Ctr  = db2Vec - mean(db2Vec);
cCOS{3}     = hCOS(db2Vec_Ctr);
cCOS_LBL{3} = 'Centered';

% Calculates the correlation
cCOS{4}     = corr(db2Vec);
cCOS_LBL{4} = 'Correlation';

% Plots the results
hFIG = figure('Position', [50 50 1250 600]);
iPlt = 0;
for iX = 1:3
    for iY = iX+1:4
        iPlt = iPlt + 1;
        subplot(2, 3, iPlt);
        plot(cCOS{iX}(:), cCOS{iY}(:), 'x')
        xlabel(cCOS_LBL{iX}); ylabel(cCOS_LBL{iY});
    end
end