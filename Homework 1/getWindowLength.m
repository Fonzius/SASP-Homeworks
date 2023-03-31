% HW1 - DAAP
% by Xinmeng Luan, Marco Bernasconi
% 31 Mar,2023

% Get window length

clc
clear
close all
% Load audio file
 [signal, Fs] = audioread('piano.wav'); 
%  [signal, Fs] = audioread('speech.wav'); 

% the Phillips-Perron (PP) test
windowlength = 2048;
M = windowlength;
win = hann(windowlength);
hop_size = floor(windowlength/4);
num_segment = ceil((length(signal)-windowlength)/hop_size + 1);
num_padding = (num_segment - 1) * hop_size + windowlength - length(signal);
signal(end+1:end+num_padding) = 0;
start_index =zeros(num_segment,1);
end_index =zeros(num_segment,1);
signal_reshape =zeros(M,num_segment);

length_tot_signal = length(signal);


for i=1:num_segment
    start_index(i) = (i-1)*hop_size + 1;
    end_index(i)= start_index(i) + windowlength - 1; 
    signal_reshape(:,i) = signal(start_index(i):end_index(i));
end
    s = win .* signal_reshape;

for i = 1:num_segment
    [h, pValue, ~, ~] = pptest(s(:,i));
    if h == 0
        fprintf('Segment %d is stationary with p-value %f\n', i, pValue);
    else
        fprintf('Segment %d is non-stationary with p-value %f\n', i, pValue);
    end
end

