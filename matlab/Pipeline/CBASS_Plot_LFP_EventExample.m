function CBASS_Plot_LFP_EventExample(db2LFP, inAnchor, db1WinSec, inSampleRate, cEVENT, cEVENT_LABEL)
%Plotting utility of the bout pipeline to plot LFP event triggered averages

% Handles optional arguments
narginchk(4, 6)
if nargin < 5, cEVENT = {}; end
if nargin < 6, cEVENT_LABEL = {}; end
if length(cEVENT) ~= length(cEVENT_LABEL)
    sprintf('Events did not all have an assigned label');
    cEVENT_LABEL = cellfun(@(x) sprintf('Event %d', x), num2cell(1:length(cEVENT)), ...
        'UniformOutput', false);
end

% Gets the number of channels and samples of the LFP
[inNChan, inNSample] = size(db2LFP);

% Checks that the event vector has the right length
if numel(inAnchor) > 1, error('in1Event must a single positive index'); end
if mod(inAnchor, 1) ~= 0 || inAnchor < 1 || inAnchor > inNSample
    error('in1Event must a single positive index')
end

% Plots exemple of the traces around bout onset
in1Ex_RelIdx    = round(inSampleRate * db1WinSec(1)): round(inSampleRate .* db1WinSec(2));
db1Time         = in1Ex_RelIdx ./ inSampleRate;

% Selection the wanted chunk of the LFP and figures out its scale
in1Sel      = in1Ex_RelIdx + inAnchor;
db2LFP_Sel  = db2LFP(:, in1Sel);
dbScale     = 0.7 * max(range(db2LFP_Sel, 2));
db1YL       = dbScale * [-.5 inNChan + 1.5];

% Plots the LFP
hold on
for iChan = 1:inNChan
    plot(db1Time, (dbScale * (inNChan - iChan + 1)) + db2LFP(iChan, in1Sel), 'k');
end
plot([0 0], db1YL, 'k--', 'LineWidth', 1), hold on
set(gca, 'YTick', dbScale * (1:inNChan), 'YTickLabel', inNChan:-1:1);
ylabel('Channel'); xlabel('Time (s)');

% Get the colormap for event coloring
db2Jet  = jet;

% Plots the event
for iEvt = 1:length(cEVENT)
    % Checks that the events are properly formatted
    if islogical(cEVENT{iEvt})
        if length(cEVENT{iEvt}) ~= inNSample
            fprintf('%s is not properly formated\r', cEVENT_LABEL{iEvt});
            continue
        end
        in1Event = find(cEVENT{iEvt}(in1Sel));
    else
        if any(mod(cEVENT{iEvt}, 1) ~= 0) || any(cEVENT{iEvt}) < 1 || any(cEVENT{iEvt}) > inNSample
            fprintf('%s is not properly formated\r', cEVENT_LABEL{iEvt});
            continue
        end
        in1Event = cEVENT{iEvt}(ismember(cEVENT{iEvt}, in1Sel)) - in1Sel(1) + 1;
    end
    
    % Sets the color
    dbColor = db2Jet(round(size(db2Jet, 1) * (iEvt/length(cEVENT))), :);
    
    % Plots the events
    if ~isempty(in1Event)
        hPLT(iEvt) = plot([1 1] * db1Time(in1Event(1)), db1YL, '--', 'Color', dbColor);
        for iPls = 2:length(in1Event)
            plot([1 1] * db1Time(in1Event(iPls)), db1YL, '--', 'Color', dbColor)
        end
    else
        hPLT(iEvt) = plot([1 1], [nan nan], '--', 'Color', dbColor);
    end
end
if ~isempty(cEVENT) , legend(hPLT, cEVENT_LABEL); end

% Enforces the limits of the y-axis
ylim(db1YL);