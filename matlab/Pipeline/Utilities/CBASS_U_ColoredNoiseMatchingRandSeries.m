function [dbXPhaseRandMat, db1Alpha, hFIG] = CBASS_U_ColoredNoiseMatchingRandSeries(dbXMat, inDim, inNSmpOut, inNSmpFit, blFig)
%DBXPHASERANDMAT = CBASS_U_ColoredNoiseMatchingRandSeries(DBXMAT, INDIM,
%INNSMPOUT, [INNSAMPFIT]) Generate a phase randomized signal fiting a p=1/(freq.^dbAlpha)
%model to the spectrum of DBXMAT and having number of sample INNSMPOUT
%along dimension INDIM. Based on CBASS_U_PhaseRandomize1D.

%Permutes the dimension so that the phase is randomized along the dimension
%of choice
in1PermDim = 1:ndims(dbXMat);
in1PermDim([2, inDim]) = [inDim, 2];
dbXMat = permute(dbXMat, in1PermDim);

%Caclulates the indices of the input spectrum
[inNChan, inNSmpIn] = size(dbXMat);
in1IdxIn    = 1:(inNSmpIn/2) + 1;
blNyq1      = mod(inNSmpIn, 2) == 0; % even FFT in

% Checks for optional arguments
if nargin < 4, inNSmpFit = inNSmpIn/2; end
if nargin < 5, blFig = false; end

%Calculates the reference indices for the interpolation
in1IdxOut   = 1:(inNSmpOut/2) + 1;
in1InterIdx = 1 + ((in1IdxOut - 1) .* (floor(inNSmpIn/2)) ./ (floor(inNSmpOut/2)));
blNyq2      = mod(inNSmpOut, 2) == 0; % even FFT out

%Computes the FFT and the amplitude
db2FFT      = fft(dbXMat, [], 2);
db2Amp      = abs(db2FFT(:, in1IdxIn));
%Fits a p=1/(freq.^dbAlpha) model
if nargout > 1, db1Alpha    = nan(inNChan, 1); end
if blFig, iFig = 0; hFIG = nan(1, ceil(inNChan/16)); end
for iChan = 1:inNChan
    % Fits the amplitude
    if blFig 
        if mod(iChan, 16) == 1, iFig = iFig + 1; hFIG(iFig) = figure('Position', [50 50 1250 600]); end
        subplot(4, 4, iChan - (16 * floor((iChan - 1)/16)));
    end
    [~, db1FitPar] = FitColoredNoiseModel(db2Amp(iChan, 2:inNSmpFit), 1:inNSmpFit - 1, blFig);
    if blFig, title(sprintf('Channel %d', iChan)); end
    if nargout > 1, db1Alpha(iChan) = db1FitPar(2); end
    db2Amp(iChan, 2:inNSmpIn/2) = exp(db1FitPar(1))./(((2:inNSmpIn/2) - 1) .^ db1FitPar(2));
end

%Interpolates the amplitude to new length and randomizes the phase
db2Amp      = interp1(in1IdxIn, db2Amp', in1InterIdx)';
if blNyq2, db2Amp = db2Amp(:, 2:end - 1); else, db2Amp = db2Amp(:, 2:end); end
db2RndAng   = rand(size(db2Amp)) * 2 * pi - pi;

%Calculates the Nyquist frequency component if needed
if blNyq2
    if blNyq1, db1FFTNyq = db2FFT(:, inNSmpIn/2 + 1);
    else, db1FFTNyq = abs(db2FFT(:, (inNSmpIn - 1)/2)) .* sin(angle(db2FFT(:, (inNSmpIn - 1)/2))); end
else
    db1FFTNyq = [];
end

%Recompose a phase randomized LFP (for the reconstituted signal to be real,
%the first element of the fft must real the rest must be conjugate
%symectric; if the number of sample is even, element NSamp/2 + 1 must be
%real)
db2FFT_Rnd  = [db2FFT(:, 1) ...
    db2Amp  .* exp(1i .* db2RndAng) ...
    db1FFTNyq ...
    conj(fliplr(db2Amp)  .* exp(1i .* fliplr(db2RndAng)))]; % Create a symetric matrix
db2FFT_Rnd  = db2FFT_Rnd .* sqrt(inNSmpOut ./ inNSmpIn); % Scales the power
dbXPhaseRandMat = ifft(db2FFT_Rnd, [], 2); % Does the IFFT

%Reshapes the LFP in the correct dimension
dbXPhaseRandMat = permute(dbXPhaseRandMat, in1PermDim);
if nargout > 1, db1Alpha = permute(db1Alpha, in1PermDim); end
if nargout > 2, if ~blFig, hFIG = []; end, end

function [db1AmpFit, db1FitPar] = FitColoredNoiseModel(db1Amp, in1Frq, blFig)
% Utility to fit a 1/(f^alpha) model. The function uses a log transform to
% solve the fit with an ordinary linear least square
if nargin < 3, blFig = false; end

db1Y = log(db1Amp);
db1X = -log(in1Frq);
db1FitPar = regress(db1Y', [ones(size(db1Y')) db1X']);
db1AmpFit = exp(db1FitPar(1))./(in1Frq .^ db1FitPar(2));

% Plots test figure if needed for diagnosic
if blFig 
    plot(in1Frq, db1Amp); hold on;
    plot(in1Frq, db1AmpFit);
    legend('data', 'fit');
    set(gca, 'yscale', 'log');
end