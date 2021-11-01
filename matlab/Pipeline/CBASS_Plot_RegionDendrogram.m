function CBASS_Plot_RegionDendrogram(db2Tree, cBRANCH)
% Utility for plotting dendrograms and highlight branches of interest

inNClu = size(db2Tree, 1) + 1;
[~, ~, in1WrdIdx] = dendrogram(db2Tree, inNClu); hold on
inNBch = length(cBRANCH);
dbYL = ylim;
hPLT = nan(1, inNBch);
for iBch = 1:inNBch
    in1Idx = find(ismember(in1WrdIdx, cBRANCH{iBch}));
    hPLT(iBch) = plot(in1Idx, dbYL(1) * ones(size(in1Idx)), 'o', 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', db2Color(round(size(db2Color, 1) * iBch/inNBch), :));
end
legend(hPLT, cellfun(@(x) sprintf('B%d', x), num2cell(1:inNBch), 'UniformOutput', false),...
    'Location', 'eastoutside')
pause(.1)