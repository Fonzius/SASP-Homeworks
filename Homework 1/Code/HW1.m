% HW1 - DAAP
% by Xinmeng Luan, Marco Bernasconi
% 31 Mar,2023

clc;
clear;
close all; 

%% 
windowLength_piano = 2048; 
windowLength_speech = 2048; 
p_piano = 10;  %  get 5 peaks, p = 10
p_speech = 46;  %  44100/1000 = 44.1

% 1 for solving Wiener-Hopf equations, 2 for steepest descent method
method = 2; 
% for Transient behaviour test of steepest descent method, 1 is doing test 
test_tran = 0;

[aPiano, fs, M, num_segment,piano_fft, tot_length, start_index, end_index ] ...
    = LPCFilter("piano.wav", windowLength_piano, p_piano, method, test_tran);
[aSpeech, ~, ~, ~,speech_fft, ~, ~, ~] ...
    = LPCFilter("speech.wav", windowLength_speech, p_speech, method, test_tran);


%% Frequency domain filter

filterApiano_freq = zeros(M,num_segment);
filterHpiano_freq = zeros(M,num_segment);
filterAspeech_freq = zeros(M,num_segment);
filterHspeech_freq = zeros(M,num_segment);
piano_H_norm = zeros(M,num_segment);
speech_H_norm = zeros(M,num_segment);
error_piano =  zeros(M,num_segment);
error_speech =  zeros(M,num_segment);
talking_freq = zeros(M,num_segment);
error_piano_time = zeros(M,num_segment);
talking_time = zeros(M,num_segment);

aPiano =aPiano';
aSpeech =aSpeech';
piano_fft = piano_fft';
speech_fft = speech_fft';

% A,H filter in freq. domain
for ss = 1:num_segment
    [filterApiano_freq(:,ss),~] = freqz(aPiano(:,ss),1,"whole",M);
    [filterHpiano_freq(:,ss),~] = freqz(1,aPiano(:,ss),"whole",M);
    [filterAspeech_freq(:,ss),~] = freqz(aSpeech(:,ss),1,"whole",M);
    [filterHspeech_freq(:,ss),~] = freqz(1,aSpeech(:,ss),"whole",M);
    piano_H_norm(:,ss) = (filterHpiano_freq(:,ss)/max(abs(filterHpiano_freq(:,ss))))*max(abs(piano_fft(:,ss)));
    speech_H_norm(:,ss) = (filterHspeech_freq(:,ss)/max(abs(filterHspeech_freq(:,ss))))*max(abs(speech_fft(:,ss)));    
    error_piano(:,ss) =  piano_fft(:,ss)./ piano_H_norm(:,ss);  
    error_speech(:,ss) =  speech_fft(:,ss)./ speech_H_norm(:,ss); 
    talking_freq(:,ss) = error_piano(:,ss) .* speech_H_norm(:,ss);
    error_piano_time(:,ss) = ifft(error_piano(:,ss));
    talking_time(:,ss) = ifft(talking_freq(:,ss));
end
talking_time = talking_time';


%% Back to time domain (OLA)
matrix_synth = zeros(num_segment,tot_length);
for i = 1:length(start_index)
    matrix_synth(i,start_index(i):end_index(i)) = talking_time(i,:);
end
output_time = sum(matrix_synth,1);
output = output_time / max(abs(output_time));
%sound(real(output),fs)

%% Audiowrite

if ~exist('GenerateSound','dir')
    mkdir('GenerateSound');
end

filename = ['GenerateSound/TI-M 20000 steps' num2str(method) '-wp' num2str(windowLength_piano) '-ws' num2str(windowLength_speech) '-pp' num2str(p_piano) '-ps' num2str(p_speech) '.wav'];
audiowrite(filename,output,fs);

