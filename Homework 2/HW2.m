clc
clear all
close all
%%
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

%% Apply STFT to both signals

window_length = 1024; % Final padded segment length
hop_size = 512; % Size of the segments of the signal

reshape_h = ceil(length(y1)/hop_size);
reshape_size = prod([hop_size,reshape_h]);

y1_zeros = [y1' zeros(1, abs(size(y1,1)-reshape_size))]';
y2_zeros = [y2' zeros(1, abs(size(y2,1)-reshape_size))]';

y1_reshape = reshape(y1_zeros,hop_size, reshape_h)';
y2_reshape = reshape(y2_zeros,hop_size, reshape_h)';

y1_padded = [y1_reshape  zeros(size(y1_reshape,1), window_length-size(y1_reshape,2))];
y2_padded = [y2_reshape  zeros(size(y2_reshape,1), window_length-size(y2_reshape,2))];

y1_stft = fft(y1_padded, window_length, 2);
y2_stft = fft(y2_padded, window_length, 2);

%% Compute feature matrix

A1 = abs(y1_stft)./(sqrt(abs(y1_stft).^2 + abs(y2_stft).^2));
A2 = abs(y2_stft)./(sqrt(abs(y1_stft).^2 + abs(y2_stft).^2));
P = (1/(2*pi)) .* angle(y2_stft./y1_stft);
phi = cat(3, A1, A2, P);

%% K-Mean Clustering

A1_unwrap = reshape(A1, numel(A1), 1);
A2_unwrap = reshape(A2, numel(A2), 1);
P_unwrap = reshape(P, numel(P), 1);

phi_unwrap = cat(2, A1_unwrap, A2_unwrap, P_unwrap);

[labels_unwrap, centroids] = kmeans(phi_unwrap, 3);

labels = reshape(labels_unwrap, size(A1'))';

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
hold on
%plot3(phi_unwrap(:,1), phi_unwrap(:,2), phi_unwrap(:,3), '.', 'Color','k');


    


%% Mask Creation

mask1 = labels;
mask1(mask1 == 2 | mask1 == 3 ) = 0;
mask1(mask1 == 1) = 1;

mask2 = labels;
mask2(mask2 == 1 | mask2 == 3 ) = 0;
mask2(mask2 == 2) = 1;

mask3 = labels;
mask3(mask3 == 1 | mask3 == 2 ) = 0;
mask3(mask3 == 3) = 1;

%% Apply Masks

f1_stft = y1_stft.*mask1;
f2_stft = y1_stft.*mask2;
f3_stft = y1_stft.*mask3;

f1_reshape = ifft(f1_stft, window_length, 2);
f2_reshape = ifft(f2_stft, window_length, 2);
f3_reshape = ifft(f3_stft, window_length, 2);

f1_reshape(:, hop_size+1:end) = [];
f2_reshape(:, hop_size+1:end) = [];
f3_reshape(:, hop_size+1:end) = [];

f1 = reshape(f1_reshape', numel(f1_reshape),1);
f2 = reshape(f2_reshape', numel(f2_reshape),1);
f3 = reshape(f3_reshape', numel(f3_reshape),1);

%% Output
%sound(real(f1), fs);
audiowrite("tryyy.wav",real(f3), fs);


%% Try to istft

y1istft = ifft(y1_stft, window_length, 2);
y1dezero = y1istft;
y1dezero(:, hop_size+1:end) = [];

y1output = reshape(y1dezero',numel(y1dezero),1);

sound(y1output, fs);


%audiowrite("tryyy.wav",f1, fs);




