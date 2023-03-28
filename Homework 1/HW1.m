clc;
clear all;
close all;

%% Legend
% p: number past samples
% n: current time
% k: k_th coefficient
% M: segment length
% a: filter coefficients

% speech: 5ms, 

%%
% [hPiano, ePiano] = LPCFilter("piano.wav");
% [hSpeech, eSpeech] = LPCFilter("speech.wav");
[aPiano] = LPCFilter("piano.wav");
[aVoice] = LPCFilter("speech.wav");
%%

aPianoNorm = aPiano./aPiano(:,1);
%% TEST WITH LPC FUNCTION
piano = audioread("piano.wav");
voice = audioread("speech.wav");
M = 220;
num_segment = ceil(length(piano)/M);

num_pad = num_segment* M -length(piano);
paddedPiano = padarray(piano,[num_pad 0],0,'post');
sPiano1 = reshape(paddedPiano,M,num_segment)';
aPiano1 = lpc(sPiano1',218);
APiano1 = zeros(size(sPiano1));
for ii = 1:size(sPiano1,1)
    APiano1(ii,:) = freqz(1, aPiano1(ii,:), M);
end

paddedVoice = padarray(voice,[num_pad 0],0,'post');
sVoice1 = reshape(paddedVoice,M,num_segment)';
aVoice1 = lpc(sVoice1',218);
AVoice1 = zeros(size(sVoice1));
for ii = 1:size(sVoice1,1)
    AVoice1(ii,:) = freqz(1, aVoice1(ii,:), M);
end

HVoice1 = 1./AVoice1;

whitenedPiano = zeros(size(sPiano1));
for ii = 1:length(whitenedPiano)
    whitenedPiano(ii,:) = fft(sPiano1(ii,:));
    whitenedPiano(ii,:) = whitenedPiano(ii,:) .* APiano1(ii,:);
end

whitenedPiano = whitenedPiano.*HVoice1;

for ii = 1:length(whitenedPiano)
    whitenedPiano(ii,:) = ifft(whitenedPiano(ii,:));
end

output = reshape(whitenedPiano, size(paddedPiano));


%%

deltaA = aPianoNorm - aPiano1;
%%
piano = audioread("piano.wav");

paddedSignal = zeros(size(A));
paddedSignal(1:length(piano)) = piano(:);

zPiano = fft(paddedSignal);
N = length(zPiano);
z = exp(2*pi*1i/N);
k = 0:N-1;
zPiano = ((1/N) * zPiano' .* z.^(-k))';

A = (1/hPiano)';

whitenedPiano = A.*paddedSignal;

pianoVoice = whitenedPiano.*hSpeech;

%pianoVoice = (N * pianoVoice' ./ z.^(-k))';
pianoVoice = ifft(piano);
N = length(pianoVoice);



sound(real(pianoVoice), 44100);

%%

wPiano = LPCSteepestDescent("piano.wav");






