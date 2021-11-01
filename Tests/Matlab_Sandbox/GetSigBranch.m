function cBRANCH = GetSigBranch(db2Tree, in1Sig)

% Creates a reference for the indices of the nodes
inNLeaves   = size(db2Tree, 1) + 1;

% Initializes a cell array grouping the clusters belinging to specific
% branches
cBRANCH = {};
% Finds branches of the dendrogram where all the leaves are significant
while ~isempty(in1Sig)
    in1Leaves       = in1Sig(1);
    inCurrentBranch = in1Sig(1);
    blAllSig        = true;
    while blAllSig
        [inNode, inLinked]  = FindLinked(db2Tree, inCurrentBranch);
        in1NewLeaves        = FindLeaves([], db2Tree, inLinked);
        if all(ismember(in1NewLeaves, in1Sig))
            in1Leaves = cat(2, in1Leaves, in1NewLeaves);
            inCurrentBranch = inNLeaves + inNode; 
        else
            blAllSig = false;
        end
    end
    cBRANCH = cat(2, cBRANCH, {in1Leaves});
    in1Sig(ismember(in1Sig, cBRANCH{end})) = [];
end

function [inNode, inLinked] = FindLinked(db2Tree, inCurrentBranch)
inNode = find(db2Tree(:, 1) == inCurrentBranch);
if ~isempty(inNode), inLinked = db2Tree(inNode, 2); return; end
inNode = find(db2Tree(:, 2) == inCurrentBranch);
inLinked = db2Tree(inNode, 1);

function in1Leaves  = FindLeaves(in1Leaves, db2Tree, inStartNode)
if inStartNode > size(db2Tree, 1) + 1
    in1Leaves  = FindLeaves(in1Leaves, db2Tree, db2Tree(inStartNode - size(db2Tree, 1) - 1, 1));
    in1Leaves  = FindLeaves(in1Leaves, db2Tree, db2Tree(inStartNode - size(db2Tree, 1) - 1, 2));
else
    in1Leaves = cat(2, in1Leaves, inStartNode);
end