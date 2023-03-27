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
[predict_piano] = LPCFilter("piano.wav");

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













