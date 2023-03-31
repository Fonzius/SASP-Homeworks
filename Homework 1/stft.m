% HW1 - DAAP
% by Xinmeng Luan, Marco Bernasconi
% 31 Mar,2023

% Get stft of the audio file
clc
clear
close all
[x, fs] = audioread('piano.wav');


%% .mp4

% STFT parameters
win_size = 1024;
overlap_pct = 50;

% Calculate STFT
[s, f, t] = spectrogram(x, win_size, round(overlap_pct/100*win_size), win_size, fs);

% Plot initial spectrogram
fig = figure();
plot_handle = plot(f, 20*log(abs(s(:, 1))));
axis([0 max(f) 0 max(20*log(abs(s(:))))]);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
xlim([0 2500])
title('Spectrogram');
grid on;

% Set up video writer object
writerObj = VideoWriter('dbpiano-spectrogram.mp4', 'MPEG-4');
writerObj.FrameRate = 30;
open(writerObj);

% Update spectrogram over time and save frames to video
for i = 2:length(t)
    set(plot_handle, 'YData', abs(s(:, i)));
    title(['Time = ', num2str(t(i)), ' s']);
    drawnow();
    pause(0.01);
    frame = getframe(fig);
    writeVideo(writerObj, frame);
end

% Close video writer object
close(writerObj);

