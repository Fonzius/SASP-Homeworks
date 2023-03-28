clc;
clear all;
close all;
%%

[piano] = audioread("piano.wav");
[voice, fs] = audioread("speech.wav");


M = 220;
for ii = 1:M:length(piano)
    a(ii,:) = weinerLPC(piano(ii:ii+M-1));
end