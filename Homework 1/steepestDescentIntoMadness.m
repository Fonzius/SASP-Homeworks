clc;
clear all;
close all;
%% 
[piano,fs] = audioread("piano.wav");

M = 20;

n_segments = ceil(length(piano)/M);
num_pad = n_segments* M -length(piano);
piano = padarray(piano,[num_pad 0],0,'post');
u = reshape(piano,M,n_segments)';

%u = reshape(piano, [M n_segments])';

p = zeros(size(u,1), size(u,2)-1);
R = zeros(M-1, M-1, n_segments);
R_max_eigen = zeros(n_segments,1);
R_min_eigen = zeros(n_segments,1);

% Compute all the mecessary stuff
for ii = 1:n_segments
    r = xcorr(u(ii,:));
    p(ii,:) = r(ceil((end+1)/2)-1:-1:1); 
    R(:,:,ii) = toeplitz(p(ii,:));
    R_eigen = eig(R(:,:,ii)); 
    R_max_eigen(ii) = max(R_eigen); 
    R_min_eigen(ii) = min(R_eigen); %Here sometimes it's negative. I don't know if it's supposed to do that
end


%%
mu = 1./R_max_eigen; %taken from the middle of the range
tau = abs(ceil(1./(2.*mu.*R_min_eigen))); %they are usually too big
w = zeros(size(u,1),size(u,2)-1);

for ii = 1:n_segments
    for jj = 1:tau(ii)
        w(ii,:) = w(ii,:) + mu(ii).*(p(ii,:)-(R(:,:,ii)*w(ii,:)')'); %formula from the slides
    end
end


%%
a_ones = ones(size(w,1),1);
a = [a_ones -w];
% A = zeros(size(a));
% for ii = 1:n_segments
%     A(ii,:) = freqz(a(ii,:),1,M);
% end

%% Test to check against

pred = zeros(n_segments, M);
for ss = 1:n_segments
    [pred(ss,:) err(ss,:)] = lpc(u(ss,:), M-1);
end

%% Difference between our results and lps funcrtion

delta = pred-a;

%% Filtering, it's wrong

H = 1./A;

u_fft = (fft(u'))';

u_wh = A.*u_fft;

u_wh_audio = (ifft(u_wh'))';

u_wh_audio_reshape = reshape(u_wh_audio, size(piano));

sound(abs(u_wh_audio_reshape), 44100);
%%
for ii = 1:n_segments
    u_wh_audio(ii,:) = filter(a(ii,:), 1, u(ii,:));
end
u_wh_audio_reshape = reshape(u_wh_audio, size(piano));
%%

audioOut =zeros(size(u));
for ss = 1:n_segments
    audioOut(ss,:) = filter([0 -a(ss,2:end)], 1, u(ss,:));
end

audioOut_reshape = reshape(audioOut',size(piano));
sound(u_wh_audio_reshape, 44100);

