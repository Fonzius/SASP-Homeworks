clc
clear all
close all

%% Data
M = 2; % Number of Microphones
N = 3; % Number of Speech Signals
d = 9e-2; %Distance between microphones
fs = 8e3; %Sampling frequency 8kHz
theta = [30, 85, -40].* (pi/180); % Angular Location of the Sources. Normal to the line between the microphones
d1 = 75e-2; %Distance of the sources from first microphone
d2 = sqrt((d-d1.*cos(theta)).^2 + d1.*sin(theta).^2); %Distance of the sources from second microphone
y1 = audioread("y1.wav");
y2 = audioread("y2.wav");

% Noise is IID (Indipendent, Identically Distributed)
sigma = 10e-3; % Std deviation of the noise of the sensor
c = 340; %Speed of sound


%% COLA STFT of the signals 

window_length = 512;
hop_size = 256;
padded_length = 1024;

% Setup for ola
window = hann(window_length)';
reshape_h = ceil(length(y1)/hop_size);
reshape_size = prod([window_length,reshape_h]);


% Pad with zeros to do the reshape
y1_zeros = [y1' zeros(1, abs(size(y1,1)-reshape_size))]';
y2_zeros = [y2' zeros(1, abs(size(y2,1)-reshape_size))]';

% Split the original signal into overlapping segments
y1_reshape = zeros(reshape_h, window_length);
y2_reshape = zeros(reshape_h, window_length);

for ii = 1:reshape_h
    index = (ii-1)*hop_size+1;
    y1_reshape(ii,:) = y1_zeros(index:index+window_length-1);
    y2_reshape(ii,:) = y2_zeros(index:index+window_length-1);
end

% Apply the window to the segments
y1_window = y1_reshape.*window;
y2_window = y2_reshape.*window;

% Pad the segments to length padded_length
y1_padded = [y1_window  zeros(size(y1_window,1), padded_length-size(y1_window,2))];
y2_padded = [y2_window  zeros(size(y2_window,1), padded_length-size(y2_window,2))];

% Take fft of all the segments to obtain the stft
y1_stft = fft(y1_padded, padded_length, 2);
y2_stft = fft(y2_padded, padded_length, 2);

% Take half of the frequencies to avoid problems with k-means
y1_stft_half = y1_stft(:, 1:padded_length/2);
y2_stft_half = y2_stft(:, 1:padded_length/2);


%% Compute feature matrix

A1 = abs(y1_stft_half)./sqrt(abs(y1_stft_half).^2 + abs(y2_stft_half).^2);
A2 = abs(y2_stft_half)./sqrt(abs(y1_stft_half).^2 + abs(y2_stft_half).^2);
P = 1/(2*pi) * angle(y2_stft_half./y1_stft_half);

%% K-Mean Clustering

% Turn the matrices into vectors and concatenate them to do the k-mean
% and then turn the labels back to the same matrix form to obtain the
% labels for each time/frequency point
A1_unwrap = reshape(A1, numel(A1), 1);
A2_unwrap = reshape(A2, numel(A2), 1);
P_unwrap = reshape(P, numel(P), 1);

phi_unwrap = cat(2, A1_unwrap, A2_unwrap, P_unwrap);

[labels_unwrap, centroids] = kmeans(phi_unwrap, 3);
labels = reshape(labels_unwrap, size(A1));

%% 3D-Plot of features

m1 = labels_unwrap;
m1(m1 == 2 | m1 == 3) = 0;
k1 = m1.* phi_unwrap;

m2 = labels_unwrap;
m2(m2 == 1 | m2 == 3) = 0;
m2(m2 == 2) = 1;
k2 = m2.* phi_unwrap;

m3 = labels_unwrap;
m3(m3 == 1 | m3 == 2) = 0;
m3(m3 == 3) = 1;
k3 = m3.* phi_unwrap;

figure()
plot3(k1(:,1), k1(:,2), k1(:,3), '.', 'Color','r');
hold on
plot3(k2(:,1), k2(:,2), k2(:,3), '.', 'Color','g');
hold on
plot3(k3(:,1), k3(:,2), k3(:,3), '.', 'Color','b');
xlabel('A1');
ylabel('A2');
zlabel('P');
legend('f1', 'f2', 'f3');


%% Mask Creation

% Turn the labels matrix into the masks for every signal
mask1 = labels;
mask1(mask1 == 2 | mask1 == 3 ) = 0;
mask1(mask1 == 1) = 1;

mask1_mirrored = [mask1, flip(mask1,2)];

mask2 = labels;
mask2(mask2 == 1 | mask2 == 3 ) = 0;
mask2(mask2 == 2) = 1;

mask2_mirrored = [mask2, flip(mask2,2)];

mask3 = labels;
mask3(mask3 == 1 | mask3 == 2 ) = 0;
mask3(mask3 == 3) = 1;

mask3_mirrored = [mask3, flip(mask3,2)];

%% Apply Masks 

% Apply the filtering with the mask
f1_stft = y1_stft.*mask1_mirrored;
f2_stft = y1_stft.*mask2_mirrored;
f3_stft = y1_stft.*mask3_mirrored;

% Take the inverse fft to obtain the segments
f1_reshape = ifft(f1_stft, padded_length, 2);
f2_reshape = ifft(f2_stft, padded_length, 2);
f3_reshape = ifft(f3_stft, padded_length, 2);

% Sum all the segments back together with overlap to obtain 
% the resulting signal
f1 = zeros(length(y1)+2*padded_length, 1);
f2 = zeros(length(y1)+2*padded_length, 1);
f3 = zeros(length(y1)+2*padded_length, 1);


for ii = 1:reshape_h
    index = ii*hop_size;
    f1(index:index+padded_length-1) = f1(index:index+padded_length-1) + f1_reshape(ii,:)';
    f2(index:index+padded_length-1) = f2(index:index+padded_length-1) + f2_reshape(ii,:)';
    f3(index:index+padded_length-1) = f3(index:index+padded_length-1) + f3_reshape(ii,:)';
end


%% Output
f1_true = audioread("DAAP_HW2_reference_files/s1.wav");
f2_true = audioread("DAAP_HW2_reference_files/s2.wav");
f3_true = audioread("DAAP_HW2_reference_files/s3.wav");

f = [f1,f2,f3];
f_true = [f1_true, f2_true, f3_true];

maxCorr = zeros(3,3);
for ii = 1:3
    for jj = 1:3
        maxCorr(ii,jj) = max(abs((xcorr(f(:,ii),f_true(:,jj)))));
    end
end
[M,MaxIndexes] = max(maxCorr,[],2);

f1_reconstructed = removeZeros(f(:,MaxIndexes(1)));
f2_reconstructed = removeZeros(f(:,MaxIndexes(2)));
f3_reconstructed = removeZeros(f(:,MaxIndexes(3)));


audiowrite("f1.wav",real(f1_reconstructed), fs);
audiowrite("f2.wav",real(f2_reconstructed), fs);
audiowrite("f3.wav",real(f3_reconstructed), fs);

%% Plots

% Log-log-amplitude spectrograms for the mixture signals (subplot).
f = linspace(0,fs/2,padded_length/2)';
t = linspace(0,length(y1)/fs, reshape_h)';


figure()
sgtitle("Log-log-amplitude spectrograms for the mixture signals");
subplot(1,2,1);
surf(t,f,20*log10(abs(y1_stft_half')), EdgeColor="none");
title("y1.wav");
xlabel("Time [s]");
ylabel("Frequency [Hz]");
view([0,90]);

subplot(1,2,2);
surf(t,f,20*log10(abs(y2_stft_half')), EdgeColor="none");
title("y2.wav");
xlabel("Time [s]");
ylabel("Frequency [Hz]");
view([0,90]);


%% Log-log-amplitude spectrograms for the true and estimated signals

f1_true = audioread("DAAP_HW2_reference_files/s1.wav");
f2_true = audioread("DAAP_HW2_reference_files/s2.wav");
f3_true = audioread("DAAP_HW2_reference_files/s3.wav");

f1_true_stft = mystft(f1_true, hop_size, window_length, padded_length);
f2_true_stft = mystft(f2_true, hop_size, window_length, padded_length);
f3_true_stft = mystft(f3_true, hop_size, window_length, padded_length);

f1_true_stft = f1_true_stft(:, 1:padded_length/2);
f2_true_stft = f2_true_stft(:, 1:padded_length/2);
f3_true_stft = f3_true_stft(:, 1:padded_length/2);

f_stft = cat(3,f1_stft,f2_stft,f3_stft);

f1_reconstructed_stft = f_stft(:,1:end/2,MaxIndexes(1));
f2_reconstructed_stft = f_stft(:,1:end/2,MaxIndexes(2));
f3_reconstructed_stft = f_stft(:,1:end/2,MaxIndexes(3));

figure()
sgtitle("Log-log-amplitude spectrograms for the true and estimated signals");

subplot(2,3,1);
surf(t,f,20*log10(abs(f1_true_stft')), EdgeColor="none");
title("True signal 1");
xlabel("Time [s]");
ylabel("Frequency [Hz]");
view([0,90]);

subplot(2,3,2);
surf(t,f,20*log10(abs(f2_true_stft')), EdgeColor="none");
title("True signal 2");
xlabel("Time [s]");
ylabel("Frequency [Hz]");
view([0,90]);

subplot(2,3,3);
surf(t,f,20*log10(abs(f3_true_stft')), EdgeColor="none");
title("True signal 3");
xlabel("Time [s]");
ylabel("Frequency [Hz]");
view([0,90]);

subplot(2,3,4);
surf(t,f,20*log10(abs(f1_reconstructed_stft')), EdgeColor="none");
title("Estimated signal 1");
xlabel("Time [s]");
ylabel("Frequency [Hz]");
view([0,90]);

subplot(2,3,5);
surf(t,f,20*log10(abs(f2_reconstructed_stft')), EdgeColor="none");
title("Estimated signal 2");
xlabel("Time [s]");
ylabel("Frequency [Hz]");
view([0,90]);

subplot(2,3,6);
surf(t,f,20*log10(abs(f3_reconstructed_stft')), EdgeColor="none");
title("Estimated signal 3");
xlabel("Time [s]");
ylabel("Frequency [Hz]");
view([0,90]);

%% Binary Masks in black and white

figure()
sgtitle("Binary Masks in black and white");

subplot(1,3,1);
surf(t,f,mask1', EdgeColor="none");
title("Mask 1");
xlabel("Time [s]");
ylabel("Frequency [Hz]");
view([0,90]);
colormap('gray');

subplot(1,3,2);
surf(t,f,mask2', EdgeColor="none");
title("Mask 2");
xlabel("Time [s]");
ylabel("Frequency [Hz]");
view([0,90]);
colormap('gray');

subplot(1,3,3);
surf(t,f,mask3', EdgeColor="none");
title("Mask 3");
xlabel("Time [s]");
ylabel("Frequency [Hz]");
view([0,90]);
colormap('gray');

%% Density Plot of each feature

figure(1)
title("A1/A2");
h1 = scatterhist(A1_unwrap,A2_unwrap, Group=labels_unwrap, Kernel='on');
xlabel("A1");
ylabel("A2");

hold on

plot(centroids(:,1), centroids(:,2), "Marker","+", "MarkerSize",10, LineWidth = 4, LineStyle='none');
legend("Speaker 1", "Speaker 2", "Spaker 3", "Centroids");



figure(2)
title("A1/P");
scatterhist(A1_unwrap,P_unwrap, Group=labels_unwrap, Kernel='on');
xlabel("A1");
ylabel("P");

hold on

plot(centroids(:,1), centroids(:,3), "Marker","+", "MarkerSize",10, LineWidth = 4, LineStyle='none');
legend("Speaker 1", "Speaker 2", "Spaker 3", "Centroids");


figure(3)
title("A2/P");
scatterhist(A2_unwrap,P_unwrap, Group=labels_unwrap, Kernel='on');
xlabel("A2");
ylabel("P");

hold on

plot(centroids(:,2), centroids(:,3), "Marker","+", "MarkerSize",10, LineWidth = 4, LineStyle='none');
legend("Speaker 1", "Speaker 2", "Spaker 3", "Centroids");

%% Density plots but made differently

nfaces = 20;
figure()
sgtitle("Density Plot for each feature");

subplot(1,3,1);
hist3([A1_unwrap A2_unwrap], 'CdataMode', 'auto', 'Nbins', [nfaces,nfaces]);
title("A1/A2");
view(2);
xlabel("A1");
ylabel("A2");

subplot(1,3,2);
hist3([A1_unwrap P_unwrap], 'CdataMode', 'auto', 'Nbins', [nfaces,nfaces]);
title("A1/P");
view(2);
xlabel("A1");
ylabel("P");

subplot(1,3,3);
hist3([A2_unwrap P_unwrap], 'CdataMode', 'auto', 'Nbins', [nfaces,nfaces]);
title("A2/AP");
view(2);
xlabel("A2");
ylabel("P");

%% TEST FEATURES

B = normalize(log10((abs(y1_stft_half)./abs(y2_stft_half))));
B = B./max(abs(B),[],"all");

%histogram(B, 5000);

R = zeros(size(y1_stft_half));
for ii = 1:size(y1_stft_half,2)
    r = xcorr(y1_stft_half(:,ii),y2_stft_half(:,ii), 'normalized');
    R(:,ii) = abs(r(ceil(end/2):end)./max(r));
end