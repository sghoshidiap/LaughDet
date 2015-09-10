function output = rhythm(y);
% if nargin < 1, y = wavread('technobeat.wav'); end
if nargin < 1, y = wavread('sgbb1.wav'); end

% Common Variables
maxSigFreq=44100;
%maxSigFreq=16000;
hannWinLength = 0.4;
hannlen = hannWinLength*2*maxSigFreq;

bandlimits=[0 200 400 800 1600 3200];
numBands = length(bandlimits);
bandRanges = zeros(12,1);
yFFTed = fft(y);
lenY = length(yFFTed)

output = zeros(lenY,numBands); %for storing the output
size(output)

i = [];
% Using band limits to find the ranges
for i = 1:numBands-1
    bandRanges(2*(i)-1) = floor(bandlimits(i)/maxSigFreq*lenY/2)+1;
    bandRanges(2*i) = floor(bandlimits(i+1)/maxSigFreq*lenY/2);
end
% The end cases of the band limits are special cases
bandRanges(2*numBands-1) = floor(bandlimits(numBands)/maxSigFreq*lenY/2)+1;
bandRanges(numBands*2) = floor(lenY/2);

% Using the frequency bands to format the output
for i = 1:numBands
    start = bandRanges(2*i-1);
    stop = bandRanges(2*i);
    output(start:stop,i) = yFFTed(start:stop);
    output(lenY+1-stop:lenY+1-start,i) = yFFTed(lenY+1-stop:lenY+1-start);
    size(output)
end
output(1,1)=0;

pause(5);

sigTime = ones(size(output));
sigFreq = ones(size(sigTime));

% For the Hann
hann = zeros(lenY,1);

% For the (relatively) small value of hannlen, it's slower to use the GPU.
for a = 1:hannlen 
    hann(a) = (cos(a*pi/hannlen/2)).^2;
end

% Going from frequency domain and then back with absolute values in the time domain
k = [];
for k = 1:numBands
    sigFreq(:,k) = fft(abs(real(ifft(output(:,k)))));
end

% Half-Hanning FFT * Signal FFT. And then back to time domain
i = [];

for i = 1:numBands
    output(:,i) = real(ifft(sigFreq(:,i).*fft(hann)));
end
