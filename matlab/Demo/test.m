% Formats the data
tmp = []
for idx = 1:100
    
    sREC                = CBASS_L1_AddPhaseRandomizedSignal(sREC);
    sTRGH_RND   = CBASS_L1_GetTrough(sREC.db2LFP_Rnd, sREC.inSampleRate, [15,30], 15);
    if blZScore, db2Data = zscore(sTRGH_RND.db2Trough);
    else, db2Data = sTRGH_RND.db2Trough; end

    % Cluster the data
    in1CluKM_Rnd = kmeans(db2Data, inNClu, 'MaxIter', inNMaxIter);

    % Get the filter for each enriched region
    db2FiltMat_Rnd = [];
    for iClu = 1:inNClu
        % Makes the filter
        in1EventIdx = sTRGH_RND.in1Index(in1CluKM_Rnd == iClu);
        db2Filter   = CBASS_U_MakeFilter(sREC.db2LFP_Rnd, sREC.inSampleRate, in1EventIdx, db1Filter_WinSec);

        % Aggregates the pulse events for the cluster of interest
        db2FiltMat_Rnd  = cat(2, db2FiltMat_Rnd, db2Filter(:));
    end

    %Calculates the correlations between random fitlers and sets the
    %threshold
    db2Corr_Rnd = corr(db2FiltMat_Rnd);
    db1Corr_Rnd   = []; 
    for iEv = 1:size(db2Corr_Rnd, 1) - 1
        db1Corr_Rnd = cat(1, db1Corr_Rnd, db2Corr_Rnd(iEv + 1:end, iEv));
    end
    dbThrsMean  = mean(db1Corr_Rnd);

    %Calculate the correlation between fitlters and sets those inferior
    %to the threshold to zero to built the adjacency matrix
    db2Adj = corr(db2FiltMat);
    db2Adj(db2Adj < max(dbThrsMean, 0)) = 0;

    in1CluSC    = CBASS_U_SpectralCluster(db2Adj);
    inNFlt      = max(in1CluSC);
    disp(["itr: " + num2str(idx) + ", inNflt: " + num2str(inNFlt)])

    tmp = [tmp, inNFlt];
end
