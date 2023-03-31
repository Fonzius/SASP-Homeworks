
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

[aPiano, ePiano, M, num_segment, outputZeros, start_index, end_index ] = LPCFilter_new("piano.wav",windowLength_piano,windowtype,p_piano);
[aSpeech, eSpeech, ~, ~, ~, ~, ~] = LPCFilter_new("speech.wav",windowLength_speech,windowtype,p_speech);


%% Time domain filter 
synth_talk_time_direct = zeros(size(ePiano));
for ss =1:size(ePiano,1) 
    synth_talk_time_direct(ss,:) = filter(1,aSpeech(ss,:), ePiano(ss,:)); %%%%right!
end
matrix_synth = zeros(size(ePiano,1),length(outputZeros));
for i = 1:length(start_index)
    matrix_synth(i,start_index(i):end_index(i)) = synth_talk_time_direct(i,:);
end
output_time_direct = sum(matrix_synth,1);
output = output_time_direct / max(abs(output_time_direct));
sound(output,44100)

%% Frequency domain filter



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


