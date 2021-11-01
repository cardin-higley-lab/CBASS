function db2PhaseRandMat = CBASS_U_PhaseRandomize2D(db2Mat)
%DBXPHASERANDMAT = CBASS_U_PhaseRandomize1D(DBXMAT, INDIM)
%Generate a phase randomized image having the same spectral amplitude as
%DBXMAT 

%Caclulates the indices
    % X
inNX            = size(db2Mat, 2);
in1Qrt1_X       = 2:1 + floor((inNX - 1)/2); % indices of the 1st quandrant in X
in1Qrt2_X       = 2 + ceil((inNX - 1)/2):inNX; % indices of the 1st quandrant in X
if mod(inNX, 2) == 0, inMd_X   = inNX/2 + 1; else, inMd_X = []; end % index of the mid point if it exist
	% Y
inNY            = size(db2Mat, 1);
in1Qrt1_Y       = 2:1 + floor((inNY - 1)/2); % indices of the 1st quandrant in Y 
if mod(inNY, 2) == 0, inMd_Y   = inNY/2 + 1; else, inMd_Y = []; end % index of the mid point if it exist

%Computes the FFT and the different segments of the phase randomized
db2FFT          = fft2(db2Mat);
db2AmpQrt1      = abs(db2FFT(in1Qrt1_Y, in1Qrt1_X));
db2AngQrt1      = rand(size(db2AmpQrt1)) * 2 * pi - pi;
db2FFTRnd_Qrt1  = db2AmpQrt1 .* exp(1i .* db2AngQrt1);
db2AmpQrt2      = abs(db2FFT(in1Qrt1_Y, in1Qrt2_X));
db2AngQrt2      = rand(size(db2AmpQrt2)) * 2 * pi - pi;
db2FFTRnd_Qrt2  = db2AmpQrt2 .* exp(1i .* db2AngQrt2); 
db2Amp0X        = abs(db2FFT(1, in1Qrt1_X));
db2Ang0X        = rand(size(db2Amp0X)) * 2 * pi - pi;
db2FFTRnd_0X    = db2Amp0X .* exp(1i .* db2Ang0X);
db2AmpMdX       = abs(db2FFT(inMd_Y, in1Qrt1_X));
db2AngMdX       = rand(size(db2AmpMdX)) * 2 * pi - pi;
db2FFTRnd_MdX   = db2AmpMdX .* exp(1i .* db2AngMdX);
db2Amp0Y        = abs(db2FFT(in1Qrt1_Y, 1));
db2Ang0Y        = rand(size(db2Amp0Y)) * 2 * pi - pi;
db2FFTRnd_0Y    = db2Amp0Y .* exp(1i .* db2Ang0Y);
db2AmpMdY       = abs(db2FFT(in1Qrt1_Y, inMd_X));
db2AngMdY       = rand(size(db2AmpMdY)) * 2 * pi - pi;
db2FFTRnd_MdY   = db2AmpMdY .* exp(1i .* db2AngMdY); 

%Reconstitutes the a phase randomized image)
db2FFT_Rnd = [db2FFT(1, 1) db2FFTRnd_0X ...
    db2FFT(1, inMd_X) conj(fliplr(db2FFTRnd_0X));... %line 1
    db2FFTRnd_0Y db2FFTRnd_Qrt1 ...
    db2FFTRnd_MdY db2FFTRnd_Qrt2; ... & line 2
    db2FFT(inMd_Y, 1) db2FFTRnd_MdX ...
    db2FFT(inMd_Y, inMd_X) conj(fliplr(db2FFTRnd_MdX)); ... % line 3
    conj(flipud(db2FFTRnd_0Y)) conj(rot90(db2FFTRnd_Qrt2, 2))...
    conj(flipud(db2FFTRnd_MdY)) conj(rot90(db2FFTRnd_Qrt1, 2))]; %line 4
% db2PhaseRandMat = ifft2(db2FFT_Rnd, 'symmetric');
db2PhaseRandMat = ifft2(db2FFT_Rnd);
    