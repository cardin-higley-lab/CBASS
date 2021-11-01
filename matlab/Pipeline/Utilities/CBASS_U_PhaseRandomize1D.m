function dbXPhaseRandMat = CBASS_U_PhaseRandomize1D(dbXMat, inDim)
%DBXPHASERANDMAT = CBASS_U_PhaseRandomize1D(DBXMAT, INDIM)
%Generate a phase randomized signal having the same spectral amplitude as
%DBXMAT along dimension INDIM. 

%Permutes the dimension so that the phase is randomized along the dimension
%of choice
in1PermDim = 1:ndims(dbXMat);
in1PermDim([2, inDim]) = [inDim, 2];
dbXMat = permute(dbXMat, in1PermDim);

%Caclulates the indices
inNSmp      = size(dbXMat, 2);
inRndIdx    = 2:1 + floor((inNSmp - 1)/2);
if mod(inNSmp, 2) == 0, inLastIdx   = inNSmp/2 + 1; else, inLastIdx = []; end

%Computes the FFT and the randomized phase
db2FFT      = fft(dbXMat, [], 2);
db2Amp      = abs(db2FFT(:, inRndIdx));
db2RndAng   = rand(size(db2Amp)) * 2 * pi - pi;

%Recompose a phase randomized LFP (for the reconstituted signal to be real,
%the first element of the fft must real the rest must be conjugate
%symectric; if the number of sample is even, element NSamp/2 + 1 must be
%real)
db2FFT_Rnd  = [db2FFT(:, 1) ...
    db2Amp  .* exp(1i .* db2RndAng) ...
    db2FFT(:, inLastIdx) ...)
    conj(fliplr(db2Amp)  .* exp(1i .* fliplr(db2RndAng)))];
dbXPhaseRandMat = ifft(db2FFT_Rnd, [], 2);

%Reshapes the LFP in the correct dimension
dbXPhaseRandMat = permute(dbXPhaseRandMat, in1PermDim);