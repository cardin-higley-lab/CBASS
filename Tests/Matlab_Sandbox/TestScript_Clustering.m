% Runs the data analysis script
cd 'D:\gamma_bouts\scripts'; % < -------------Edit 
run LoadDataScript_BetaGamma_HilbertFormatting.m
close all
%% Initializes clustering parameters
% Sets clustering parameters
rng(1949);  % initializes the rando number generator
in1NClu     = [10, 20, 40, 80]; %number of clusters
% in1NClu     = 40;
inNMaxIter  = 10000; %max number of iterations
iBnd        = 2;

% Formats the data
db2Data = zscore(cMOTIF{iBnd});

% Sets verbose
blVerbose = true;

% Sets version of the algorithm
inVersion = 2;
dbSigThrs = 10.^-4;

% Initializes the figure
hFIG = figure('Position', [50, 50, 1250, 600]);
db2Color = jet;
%% Loops throught cluster size
for iNCl = 1:length(in1NClu)
    % Sets cluster number
    inNClu = in1NClu(iNCl);
    
    % Print informations on the screen
    if blVerbose
        fprintf('\r------ %d CLUSTERS ---------------------------\r', inNClu)
    end
    
    % Cluster the data
    tic; [in1CluKM, db2Ctr] = kmeans(db2Data, inNClu, 'MaxIter', inNMaxIter); toc
    
    %%% Checks if clusters are enriched for running the data
    % Computes the indices of the bouts where the mouse was running
    bl1Run_Bout     = bl1Run(cTROUGH_IDX{iBnd});
    
    % Computes the overall ratio of running points in the data
    dbRate_All      = mean(bl1Run_Bout);
    
    % Initialize a cluster structure containing stats for enrichment
    sCLU = struct('inNObs', cell(1, inNClu), ... 
        'dbRate', cell(1, inNClu), ...
        'dbRate_Dev', cell(1, inNClu), ...
        'dbPVal', cell(1, inNClu));
    
    % Loops over clusters
    for iClu = 1:inNClu
        bl1Clu  = in1CluKM == iClu; % indices of the cluster
        inNObs  = sum(bl1Clu); % number of points in the cluster
        sCLU(iClu).inNObs       = inNObs;
        
        % Calculate the ratio of running points for the cluster and the p-value
        % of that ratio under H0 i.e. the hypothesis that the ratio is the same
        % that in the overall population. The p-value is calculated using a
        % binomial distribution.
        dbRate                  = mean(bl1Run_Bout(bl1Clu)); % Calculates the ratio
        dbRate_Dev              = dbRate - dbRate_All;
        sCLU(iClu).dbRate       = dbRate; % Stores the ratio
        sCLU(iClu).dbRate_Dev   = dbRate_Dev; % Stores the ratio
        sCLU(iClu).dbPVal       = 2 * (1 - normcdf(abs(dbRate_Dev)./...
            sqrt(dbRate_All * (1 - dbRate_All) ./ inNObs)));
    end
    
    % Calculates the FDR corrected p-values (%WARNING: requires fdr_bh from the
    % matlab file exchange)
    cSIG = num2cell(fdr_bh([sCLU.dbPVal]));
    [sCLU.blSig_FDR] = cSIG{:};
    
    %%% Print the results
    % Sort the cluster by ratio
    if blVerbose
        [~, in1Sort] = sort([sCLU.dbRate_Dev]);
        
        % Initialize a significance statement array
        cSTATEMENT = {'NON SIG' 'SIG'};
        
        % Loops throught clusters
        for iClu = in1Sort
            fprintf('Cluster %2d (N = %d):\tRate: %.4f\tDev: %.3f\t(p = %.4f)\t%s\r', ...
                iClu, sCLU(iClu).inNObs, sCLU(iClu).dbRate, sCLU(iClu).dbRate_Dev, sCLU(iClu).dbPVal, ...
                cSTATEMENT{sCLU(iClu).blSig_FDR + 1});
        end
    end
    
    %%% Links the centroid by distances and identify enriched branches
    % Creates a tree of the data
    db2Tree     = linkage(db2Ctr, 'ward');
    
    % Gets the branches of the dendrogram that are significant
    if inVersion == 1    
        % Creates a pool of the positives indices and sorts them by ratio
        in1Sig          = find([sCLU.dbRate_Dev] > 0 & [sCLU.blSig_FDR]);
        [~, in1SortIdx] = sort([sCLU(in1Sig).dbRate], 'ascend' );
        in1Sig          = in1Sig(in1SortIdx);
        cBRANCH = GetSigBranch(db2Tree, in1Sig);
    elseif inVersion == 2
        in1NObs = [sCLU.inNObs]; db1Rate = [sCLU.dbRate]; 
        [cBRANCH, in1NObs_Bch, db1Rate_Bch, db1PVal_Bch] = ...
            GetSigBranch_v2(db2Tree, in1NObs, db1Rate, dbSigThrs);
    end
    % Finds the number of branches
    inNBch  = length(cBRANCH);
    
    % Create a figure of the result
    subplot(2, 2, iNCl)
    [~, ~, in1WrdIdx] = dendrogram(db2Tree, inNClu); hold on
    dbYL = ylim;
    hPLT = nan(1, inNBch);
    for iBch = 1:inNBch
        in1Idx = find(ismember(in1WrdIdx, cBRANCH{iBch}));
        hPLT(iBch) = plot(in1Idx, dbYL(1) * ones(size(in1Idx)), 'o', 'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', db2Color(round(size(db2Color, 1) * iBch/inNBch), :));
    end
    legend(hPLT, cellfun(@(x) sprintf('B%d', x), num2cell(1:inNBch), 'UniformOutput', false),...
        'Location', 'eastoutside')
    title(sprintf('%d Clusters', inNClu))
    pause(.1)
    
    % Saves the output vector
    dbTemp = zeros(size(in1CluKM));
    for iBch = 1:inNBch
        dbTemp(ismember(in1CluKM, cBRANCH{iBch})) = iBch;
    end
    eval(sprintf('in1BranchIdx_%dClusters = dbTemp;', in1NClu(iNCl)));
    
end
%% Saves the figure
if inVersion == 1
    saveas(hFIG, 'TestClustering_1.png');
    save('ClusterBranchIndices_v1_15', 'in1BranchIdx*');
elseif inVersion == 2
    saveas(hFIG, 'TestClustering_2.png');
    save('ClusterBranchIndices_v2_15', 'in1BranchIdx*');
end

