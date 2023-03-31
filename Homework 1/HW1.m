
clc;
clear all;
close all;

%%%%
% Q1: check windowing and adding 
% Q2: check steepest 
%%%%
%% Legend
% p: number past samples
% n: current time
% k: k_th coefficient
% M: segment length
% a: filter coefficients

% speech: 5ms, 

%%
windowLength_piano = 2048; %2000
windowLength_speech = 2048; %220
p_piano = 10; % test, get 5 peaks, so p = 10
p_speech = 46; % 44100/1000 = 44.1
windowtype = "hann";

[aPiano, ePianotime, ePianofreq, M, num_segment,piano_fft, outputZeros, start_index, end_index ] = LPCFilter_new("piano.wav",windowLength_piano,windowtype,p_piano);
[aSpeech, eSpeechtime, eSpeechfreq, ~, ~,speech_fft, ~, ~, ~] = LPCFilter_new("speech.wav",windowLength_speech,windowtype,p_speech);


%% Time domain filter 
synth_talk_time_direct = zeros(size(ePianotime));
for ss =1:size(ePianotime,1) 
    synth_talk_time_direct(ss,:) = filter(1,aSpeech(ss,:), ePianotime(ss,:)); %%%%right!
end
matrix_synth = zeros(size(ePianotime,1),length(outputZeros));
for i = 1:length(start_index)
    matrix_synth(i,start_index(i):end_index(i)) = synth_talk_time_direct(i,:);
end
output_time_direct = sum(matrix_synth,1);
output = output_time_direct / max(abs(output_time_direct));
 sound(output,44100)


%% Frequency test new

Hpiano = piano_fft ./ePianofreq;
Hspeech = speech_fft ./eSpeechfreq;

talking_seg_freq = ePianofreq.*Hspeech;
talking_seg_time = ifft(talking_seg_freq')';



%% Frequency domain filter

filterApiano_freq = zeros(num_segment,M);
filterHpiano_freq = zeros(num_segment,M);
filterAspeech_freq = zeros(num_segment,M);
filterHspeech_freq = zeros(num_segment,M);

% errorApiano_freq = zeros(num_segment,M);
% errorHpiano_freq = zeros(num_segment,M);
% errorAspeech_freq = zeros(num_segment,M);
% errorHspeech_freq = zeros(num_segment,M);


for ss = 1:num_segment
    [filterApiano_freq(ss,:),~] = freqz(aPiano(ss,:),1,"whole",M);
    [filterHpiano_freq(ss,:),~] = freqz(1,aPiano(ss,:),"whole",M);
    [filterAspeech_freq(ss,:),~] = freqz(aSpeech(ss,:),1,"whole",M);
    [filterHspeech_freq(ss,:),~] = freqz(1,aSpeech(ss,:),"whole",M);
end


aPiano_fft = zeros(num_segment,M);
fft1_piano = zeros(num_segment,M);
ones1 = ones(size(aPiano));
for ss = 1:num_segment
    aPiano_fft(ss,:) = fft(aPiano(ss,:),M);
    fft1_piano(ss,:) = fft(ones1(ss,:),M);  
end

Afilterpianofreq = aPiano_fft ./ fft1_piano;
Afilterpianotime = ifft(Afilterpianofreq);


error_piano_freq = piano_fft' .* filterApiano_freq;
error_speech_freq = speech_fft' .* filterAspeech_freq;

talking_freq_seg = error_piano_freq .* filterHspeech_freq;
talking_time_seg = ifft(talking_freq_seg);

% COLA
matrix_synth = zeros(size(ePiano,1),length(outputZeros));
for i = 1:length(start_index)
    matrix_synth(i,start_index(i):end_index(i)) = synth_talk_time_direct(i,:);
end
output_time_direct = sum(matrix_synth,1);
output = output_time_direct / max(abs(output_time_direct));
sound(output,44100)




% for ss = 1:length()
% [A(:,nn),~] = freqz(aPiano(:,nn),1,"whole",M);

% talking_freq = hSpeech .* ePiano;
% talking_time = ifft(talking_freq);
% 
% talking_time_cola = adding(talking_time,0.5,windowLength_piano);


% sound(abs(talking_time_cola),44100)
% clear sound

%% Audiowrite
% windowLength_piano = 2000; %2000
% windowLength_speech = 2000; %220
% p_piano = 8; % test, get 5 peaks, so p = 10
% p_speech = 50; % 44100/1000 = 44.1

if ~exist('GenerateSound','dir')
    mkdir('GenerateSound');
end

filename = ['GenerateSound/TI-wp' num2str(windowLength_piano) '-ws' num2str(windowLength_speech) '-pp' num2str(p_piano) '-ps' num2str(p_speech) '.wav'];
audiowrite(filename,output,44100);


