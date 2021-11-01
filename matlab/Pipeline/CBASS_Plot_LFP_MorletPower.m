function CBASS_Plot_LFP_MorletPower(db2LFP, inSampleRate, blPlotChan, in1State, cSTATE_LABEL)
%Plotting utility of the bout pipeline to plot LFP power

% Checks for number of arguments
narginchk(2, 5)
if nargin < 3, blPlotChan = false; end
if nargin < 4, in1State = true(size(db2LFP, 2)); end
if nargin < 5, cSTATE_LABEL = cellfun(@num2str, num2cell(unique(in1State)), 'UniformOutput', false); end

% Gets the number of channels
inNChan     = size(db2LFP, 1);

% Sets the x axis
in1FrqLbl   = 2:2:120;

% Determines states
in1U_State  = unique(in1State);
inNSt       = length(in1U_State);
cPOWER      = cell(1, inNSt);

% Aggregates the power for the LFP
for iChan = 1:inNChan
    db2Power = abs(CBASS_U_MorletTransform(db2LFP(iChan, :), inSampleRate));
    for iSt = 1:inNSt
        bl1Sel = in1State == in1U_State(iSt);
        cPOWER{iSt} = cat(1, cPOWER{iSt}, mean(db2Power(:, bl1Sel), 2)'); hold on
    end
end

% Sets the colors for the plot
if inNSt > 1
    db2Jet = jet;
    in1Sel = round(0.2 * size(db2Jet, 1)): round(0.8 * size(db2Jet, 1));
    db2Jet = db2Jet(in1Sel, :);
    in1Idx = 1 + round(linspace(0, 1, inNSt) * (size(db2Jet, 1) - 1));
    db2Color = db2Jet(in1Idx, :);
else
    db2Color = [0 0 0];
end

% Plots the mean spectrum
hPLT = nan(1, inNSt);
for iSt = 1:inNSt
    if blPlotChan, plot(in1FrqLbl, 10 * log10(cPOWER{iSt}'), 'Color', db2Color(iSt, :)./max(db2Color(iSt, :))); end
    hPLT(iSt) = plot(in1FrqLbl, 10 * log10(mean(cPOWER{iSt}, 1)), 'Color', db2Color(iSt, :), 'LineWidth', 2);
end
xlabel('Frequency (Hz)'); ylabel('Power (dB)');
if inNSt > 1; legend(hPLT, cSTATE_LABEL); end
