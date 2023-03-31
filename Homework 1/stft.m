%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%/////////Get stft of the audio file/////////%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all


  [x, fs] = audioread('piano.wav');
%  [x, fs] = audioread('speech.wav');

%% Method 1: waterfall 
% win_size = 1024;
% overlap_pct = 50;
% [s, f, t] = spectrogram(x, win_size, round(overlap_pct/100*win_size), win_size, fs);
% 

% imagesc(t, f, 20*log10(abs(s)));
% axis xy;
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');


%% Method 2 -correct

% % STFT parameters
% win_size = 1024;
% overlap_pct = 50;
% 
% % Calculate STFT
% [s, f, t] = spectrogram(x, win_size, round(overlap_pct/100*win_size), win_size, fs);
% 
% % Plot initial spectrogram
% fig = figure();
% plot_handle = plot(f, abs(s(:, 1)));
% axis([0 max(f) 0 max(abs(s(:)))]);
% xlabel('Frequency (Hz)');
% ylabel('Amplitude');
% xlim([0 5000])
% title('Spectrogram');
% grid on;
% 
% % Update spectrogram over time
% for i = 2:length(t)
%     set(plot_handle, 'YData', abs(s(:, i)));
%     title(['Time = ', num2str(t(i)), ' s']);
%     drawnow();
%     pause(0.01);
% end

%% Method 3 - save video
% 
% % Load audio file
% 
% % STFT parameters
% win_size = 1024;
% overlap_pct = 50;
% 
% % Calculate STFT
% [s, f, t] = spectrogram(x, win_size, round(overlap_pct/100*win_size), win_size, fs);
% 
% % Plot initial spectrogram
% fig = figure();
% plot_handle = plot(f, abs(s(:, 1)));
% axis([0 max(f) 0 max(abs(s(:)))]);
% xlim([0 5000])
% xlabel('Frequency (Hz)');
% ylabel('Amplitude');
% title('Spectrogram');
% grid on;
% 
% % Set up video writer object
% writerObj = VideoWriter('piano.avi');
% writerObj.FrameRate = 30;
% open(writerObj);
% 
% % Update spectrogram over time and save frames to video
% for i = 2:length(t)
%     set(plot_handle, 'YData', abs(s(:, i)));
%     title(['Time = ', num2str(t(i)), ' s']);
%     drawnow();
%     pause(0.01);
%     frame = getframe(fig);
%     writeVideo(writerObj, frame);
% end
% 
% % Close video writer object
% close(writerObj);


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

