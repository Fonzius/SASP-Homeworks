%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%/////////Get window length/////////%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all
% Load audio file
 [x, Fs] = audioread('piano.wav'); % use 2000 --> test result is staionary
% [x, Fs] = audioread('speech.wav'); % use Mezza said 5ms --> 220

%% Method 1 - test mean and variance

% % Define window size and overlap size
% win_size = 1024;
% overlap = 512;
% 
% % Create Hamming window
% ham_win = hamming(win_size);
% 
% % Divide signal into overlapping segments
% segments = buffer(x, win_size, overlap);
% 
% % Apply Hamming window to each segment
% segments_win = bsxfun(@times, segments, ham_win);
% 
% % Compute mean and variance of each segment
% seg_mean = mean(segments);
% seg_var = var(segments);
% seg_win_mean = mean(segments_win);
% seg_win_var = var(segments_win);
% 
% % Verify that mean and variance of each segment do not change
% max_diff_mean = max(abs(seg_mean - seg_win_mean));
% max_diff_var = max(abs(seg_var - seg_win_var));
% 
% if max_diff_mean < eps && max_diff_var < eps
%     disp('Mean and variance of each segment do not change.')
% end

%% Method 2 - the Augmented Dickey-Fuller (ADF) test

% % Define window size and overlap size
% win_size = 1024;
% overlap = win_size/2;
% 
% % Create Hamming window
% ham_win = hamming(win_size);
% 
% % Divide signal into overlapping segments
% segments = buffer(x, win_size, overlap);
% 
% % Apply Hamming window to each segment
% segments_win = bsxfun(@times, segments, ham_win);
% 
% % Test each segment for stationarity using ADF test
% num_segments = size(segments_win, 2);
% for i = 1:num_segments
%     [h, pValue, ~, ~] = adftest(segments_win(:,i), 'lags', floor(length(segments_win(:,i))/2));
%     if h == 0
%         fprintf('Segment %d is stationary with p-value %f\n', i, pValue);
%     else
%         fprintf('Segment %d is non-stationary with p-value %f\n', i, pValue);
%     end
% end

%% Method 3 - the Kwiatkowski-Phillips-Schmidt-Shin (KPSS) test

% % Define window size and overlap size
% win_size = 8; %1024
% overlap = win_size/2;
% 
% % Create Hamming window
% ham_win = hamming(win_size);
% 
% % Divide signal into overlapping segments
% segments = buffer(x, win_size, overlap);
% 
% % Apply Hamming window to each segment
% segments_win = bsxfun(@times, segments, ham_win);
% 
% % Test each segment for stationarity using KPSS test
% num_segments = size(segments_win, 2);
% for i = 1:num_segments
%     [h, pValue, ~, ~] = kpsstest(segments_win(:,i));
%     if h == 0
%         fprintf('Segment %d is stationary with p-value %f\n', i, pValue);
%     else
%         fprintf('Segment %d is non-stationary with p-value %f\n', i, pValue);
%     end
% end


%% method 4 : the Phillips-Perron (PP) test

% Define window size and overlap size
win_size = 2000; %1024
overlap = win_size/2;

% Create Hamming window
ham_win = hamming(win_size);

% Divide signal into overlapping segments
segments = buffer(x, win_size, overlap);

% Apply Hamming window to each segment
segments_win = bsxfun(@times, segments, ham_win);

% Test each segment for stationarity using Phillips-Perron (PP) test
num_segments = size(segments_win, 2);
for i = 1:num_segments
    [h, pValue, ~, ~] = pptest(segments_win(:,i));
    if h == 0
        fprintf('Segment %d is stationary with p-value %f\n', i, pValue);
    else
        fprintf('Segment %d is non-stationary with p-value %f\n', i, pValue);
    end
end

