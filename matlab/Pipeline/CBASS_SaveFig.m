function CBASS_SaveFig(chPath, hFIG, cFIGNAME, chOption)
% Bout pipeline utility to save figures. Save figures in the handle vector
% hFIG in the folder provided by chPath under the names provided in the
% cell array cFIGNAME; By default a .eps, .png and .pdf files are saved.
% This can be controled by the option string chOption. To save only .png
% and .pdf. Set chOption to 'pngpdf'. The routine requires exportfig from
% the file exchange.

% Checks the number of arguments
narginchk(3, 4);

% Sets the default option
chDefault = '';

% Sets the option to default if not provided
if nargin < 4, chOption = chDefault; end

% Turns off warnings
warning('off');

% Checks that their is a name for each figure in hFIG
if length(hFIG) ~= length(cFIGNAME)
    error('Provide a name for each figure');
end

% Checks the option
blEPS = contains(chOption, 'eps', 'IgnoreCase', true);
blPNG = contains(chOption, 'png', 'IgnoreCase', true);
blPDF = contains(chOption, 'pdf', 'IgnoreCase', true);
if ~any([blEPS, blPNG, blPDF]); [blEPS, blPNG, blPDF] = deal(true); end

% Saves the figures
for iFig = 1:length(hFIG)
    if blEPS
        try
            exportfig(hFIG(iFig) ,fullfile(chPath, cFIGNAME{iFig}), ...
                'Color', 'rgb', 'Renderer', 'painters');
        catch ME
            fprintf('.eps could not be saved\r');
            getReport(ME);
        end
    end
    
    if blPNG
        try
            saveas(hFIG(iFig), fullfile(chPath, [cFIGNAME{iFig} '.png']))
        catch ME
            fprintf('.png could not be saved\r');
            getReport(ME);
        end
    end
    
    if blPDF
        try
            set(hFIG(iFig), 'Units', 'centimeters')
            db1Pos = get(hFIG(iFig), 'Position');
            set(hFIG(iFig), 'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', 'PaperSize', [db1Pos(3) db1Pos(4)])
            print(hFIG(iFig), fullfile(chPath, cFIGNAME{iFig}), '-dpdf')
        catch ME
            fprintf('.pdf could not be saved\r');
            getReport(ME); 
        end
    end
end

% Turns warnings back on
warning('on');