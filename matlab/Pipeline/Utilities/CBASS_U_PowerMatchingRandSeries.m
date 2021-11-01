function dbXPhaseRandMat = CBASS_U_PowerMatchingRandSeries(dbXMat, inDim, inNSmpOut)
%DBXPHASERANDMAT = CBASS_U_PowerMatchingRandSeries(DBXMAT, INDIM, INNSMPOUT)
%Generate a phase randomized signal having the same spectral amplitude as
%DBXMAT and having number of sample INNSMPOUT along dimension INDIM.
%Based on CBASS_U_PhaseRandomize1D.

%Permutes the dimension so that the phase is randomized along the dimension
%of choice
in1PermDim = 1:ndims(dbXMat);
in1PermDim([2, inDim]) = [inDim, 2];
dbXMat = permute(dbXMat, in1PermDim);

%Caclulates the indices of the input spectrum
inNSmpIn    = size(dbXMat, 2);
in1IdxIn    = 1:(inNSmpIn/2) + 1;
blNyq1      = mod(inNSmpIn, 2) == 0; % even FFT in

%Calculates the reference indices for the interpolation
in1IdxOut   = 1:(inNSmpOut/2) + 1;
in1InterIdx = 1 + ((in1IdxOut - 1) .* (floor(inNSmpIn/2)) ./ (floor(inNSmpOut/2)));
blNyq2      = mod(inNSmpOut, 2) == 0; % even FFT out

%Computes the FFT and the amplitude
db2FFT      = fft(dbXMat, [], 2);
db2Amp      = interp1(in1IdxIn, abs(db2FFT(:, in1IdxIn))', in1InterIdx)';
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