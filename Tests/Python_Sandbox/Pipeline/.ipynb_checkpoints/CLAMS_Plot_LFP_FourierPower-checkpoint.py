import numpy as np

def Plot_LFP_FourierPower(db2LFP, inSampleRate, blPlotChan, in1State, sOPTION, cSTATE_LABEL):
    #Plotting utility of the bout pipeline to plot LFP power

    # Checks for number of arguments
#     narginchk(2, 5)
#     if nargin < 3, blPlotChan = false; end
#     if nargin < 4, in1State = true(1, size(db2LFP, 2)); end
    verbose = sOPTION.blVerbose
#     cSTATE_LABEL = cellfun(@num2str, num2cell(unique(in1State)), 'UniformOutput', false); end
    cSTATE_LABEL = str(np.unique(in1State))
    if verbose: print('cSTATE_LABEL: ',cSTATE_LABEL)

    # Gets the number of channels
    inNChan     = db2LFP.shape[0]

    # Sets the x axis
    dbWinLenSec = .5 #Default of CBASS_U_FourierPower
    inTopFreq   = 120 #Maximum frequency to plot in Hz
    in1FrqIdx   = np.arange(inTopFreq*dbWinLenSec)
    if verbose: 
        print('in1FrqIdx.shape: ', in1FrqIdx.shape)
        print('in1FrqIdx: ', in1FrqIdx)
    in1FrqLbl   = np.arange(0,inTopFreq,1/dbWinLenSec)
    if verbose: 
        print('in1FrqLbl.shape: ', in1FrqLbl.shape)
        print('in1FrqLbl: ', in1FrqLbl)

#     # Aggregates the power for the LFP
#     db3Pwr = [];
#     for iChan in range(inNChan):
#         [db2Power, in1CntrIdx] = CBASS_U_FourierPower(db2LFP(iChan, :), inSampleRate);
#         db3Pwr = cat(3, db3Pwr, db2Power);
#     end

#     % Determines states
#     in1BinState = in1State(in1CntrIdx);
#     in1U_State  = unique(in1BinState);
#     inNSt       = length(in1U_State);

#     % Sets the colors for the plot
#     if inNSt > 1
#         db2Jet = jet;
#         in1Sel = round(0.2 * size(db2Jet, 1)): round(0.8 * size(db2Jet, 1));
#         db2Jet = db2Jet(in1Sel, :);
#         in1Idx = 1 + round(linspace(0, 1, inNSt) * (size(db2Jet, 1) - 1));
#         db2Color = db2Jet(in1Idx, :);
#     else
#         db2Color = [0 0 0];
#     end

#     % Plots the mean spectrum
#     hPLT = nan(1, inNSt);
#     for iSt = 1:inNSt
#         bl1Sel = in1BinState == in1U_State(iSt);
#         db2Pwr = permute(mean(db3Pwr(in1FrqIdx, bl1Sel, :), 2), [1 3 2]); hold on
#         if blPlotChan, plot(in1FrqLbl, 10 * log10(db2Pwr), 'Color', .5 * db2Color(iSt, :) + .5); end
#         hPLT(iSt) = plot(in1FrqLbl, 10 * log10(mean(db2Pwr, 2)), 'Color', db2Color(iSt, :), 'LineWidth', 2);
#     end
#     xlabel('Frequency (Hz)'); ylabel('Power (dB)');
#     if inNSt > 1; legend(hPLT, cSTATE_LABEL); end
