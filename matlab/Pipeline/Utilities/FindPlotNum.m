function [inNRow, inNCol] = FindPlotNum(inNPlt)
% Little utility to organize suplots in a figure as a function of the
% number of plots we have.
inNCol = ceil(sqrt(inNPlt));
inNRow = round(sqrt(inNPlt));
