clc;
clear all;
close all;
%%

[piano, pianoFs] = audioread("piano.wav");
[speech, speechFs] = audioread("speech.wav");

if pianoFs == speechFs
    fs = pianoFs;
end

