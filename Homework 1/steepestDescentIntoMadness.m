clc;
clear all;
close all;
%%
[piano,fs] = audioread("piano.wav");

M = 25;
n_segments = length(piano)/M;

u = reshape(piano, [M n_segments])';

p = zeros(size(u,1), size(u,2)-1);
R = zeros(M-1, M-1, n_segments);
R_max_eigen = zeros(n_segments,1);
R_min_eigen = zeros(n_segments,1);

for ii = 1:n_segments
    r = xcorr(u(ii,:),"normalized");
    p(ii,:) = r(ceil((end+1)/2)-1:-1:1);
    R(:,:,ii) = toeplitz(p(ii,:));
    R_eigen = eig(R(:,:,ii));
    R_max_eigen(ii) = max(R_eigen);
    R_min_eigen(ii) = min(R_eigen);
end

mu = 1./R_max_eigen;
tau = ceil(1./(2.*mu.*R_min_eigen));
w = zeros(size(u,1),size(u,2)-1);

for ii = 1:n_segments
    for jj = 1:tau(ii)
        w(ii,:) = w(ii,:) + mu(ii).*(p(ii,:)-(R(:,:,ii)*w(ii,:)')');
    end
end

%%
a_ones = ones(size(w,1),1);
a = [a_ones -w];
A = zeros(size(a));
for ii = 1:n_segments
    A(ii,:) = freqz(a(ii,:),1,M);
end
%%

H = 1./A;

u_fft = (fft(u'))';

u_wh = A.*u_fft;

u_wh_audio = (ifft(u_wh'))';

u_wh_audio_reshape = reshape(u_wh_audio, size(piano));

sound(abs(u_wh_audio_reshape), 44100);