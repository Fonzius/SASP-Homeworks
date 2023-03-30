function[H_reshape, H_norm_reshape error_freq] = LPCFilter (audioFile)

% Import the files
[signal, fs] = audioread(audioFile);

%%%% add windowing to the signal
% segment == window 
% 5 ms for SPEECH! is taken from lesson as example segment length
% for harmonic and noise signal, probably we could have longer windows (
% more segments in each window)


%%%% staionary assumption, recording 20230317, 31:57 --> the mean and the
%%%% variance don't change with time. 
%%%% we need to compute the mean and variance to make sure they don't
%%%% change in each segment.
M = floor(5e-3*fs); % How many samples in each segment


%%%% recording 20230317, 39:46 --> error signal (in time) is longer than
%%%% each windowed signal, due to the convolution
%%%% p is the length of the filter
%%%% e.g. p=4 --> 2 resonant peaks
%%%%      p=8 --> 4 resonant peaks...
%%%% cuz poles are complicated conjugate

%%%% using lower order p is able to get the peak in spectrogram (for
%%%% piano), like fs/1000 <= p <= fs/1000+4 
%%%% the code works well to get the peaks instead of the valleys


%%%% get pitch from piano (should be only containing fundamental frequency)
%%%% and get spectrogram from audio

%%%% find a to minimize the short-time mean-squared error

%%%%//// Method by Xinmeng
num_segment = ceil(length(signal)/M);
num_pad = num_segment* M -length(signal);
paddedSignal = padarray(signal,[num_pad 0],0,'post');
s = reshape(paddedSignal,M,num_segment)';
s_fft = fft(paddedSignal);

s_fft_seg = fft(s')';
%%%%%////Method by Marco 
%  index = 1:M:length(signal); %resize the piano signal
% s = zeros(length(index), M); %piano signal split in segments
% 
% % Create a new signal by appending zeros at the end s.t. length(paddedSignal) is
% % multiple of M
% paddedSignal = zeros(numel(s),1);
% paddedSignal(1:length(signal)) = signal(:);
% 
% % split the signal in segments of length M
% for ii = 1:length(s)
%     s(ii,:) = paddedSignal(index(ii) : index(ii)+M-1);
% end


% 
% % Create the r vector of the autocorrelations with sample lag 1:M
% r = zeros(length(s), size(s,2)); % n , p
% for nn = 1:length(s)
%     for pp = 1:M
%         r(nn,pp) = sum(s(nn,1:M-pp).*s(nn,1+pp:M));
%     end
% end
% 
% 
% % Create autocorrelation vector with sample lags 0:M-1
% acorrVector = zeros(length(s), size(s,2));
% for nn = 1:length(s)
%     for pp = 0:M-1
%         acorrVector(nn,pp+1) = sum(s(nn,1:M-pp).*s(nn,1+pp:M));
%     end
% end




%%%%%//// auto-correlation matrix by xinmeng

auto_para_cal = zeros(M,M);
r_auto_correlation = zeros(M);
for ss = 1:num_segment
    sample_segment_current = s(ss,:)';
    % get calculation parameter matrix of auto-correlation
    for mm = 1:M        
        auto_para_cal(mm, 1:M-mm+1) = sample_segment_current(mm:end)'; % fill the matrix 
    end
    auto_correlation_current = auto_para_cal * sample_segment_current;
    r_auto_correlation(:,ss) = auto_correlation_current';
end

% test for normalization
r_auto_correlation_norm = r_auto_correlation ./ max(r_auto_correlation, [], 1);
% r_auto_correlation = r_auto_correlation';
% r_auto_correlation_norm =r_auto_correlation_norm';

r = r_auto_correlation(2:end,:);
R = zeros(M-1,M-1,num_segment);
for ss= 1:num_segment
    for ii = 1:M-1
        for jj = 1:M-1
            R(ii,jj,ss) = r_auto_correlation(abs(ii-jj)+1,ss);
        end 
    end 
end
a = zeros(M-1,num_segment);
R_inverse = zeros(M-1,M-1,num_segment);

%%%%%%%%%%% VERY SLOW
for ss= 1:num_segment
    R_inverse(:,:,ss) = inv(R(:,:,ss));
    a(:,ss) = R_inverse(:,:,ss) * r(:,ss);
end

r = r';
a = a';
a_exp1 = ones(size(a,1),1);
a_exp =[a_exp1 -1.*a];

%%% get filter

%%%%%%%%%%%% method 1 - SHIT！ 
H = zeros(size(s));
% A = zeros(size(s));
% H_1 = zeros(size(s));
for ss = 1:num_segment
    H_index = freqz(1, a_exp(ss,:),"whole",M);
    H(ss,:) = H_index';

%     A_index = freqz(a_exp(ss,:),1,"whole",M);
%     A(ss,:) = A_index';

    %%%test
%     pred_index_1 = freqz(1, a_exp(ss,:),M);
%     H_1(ss,:) = pred_index_1';
    %%%//////
%     H_current = fftshift(pred_index);
%     H(ss,:) = H_current';
end

A = 1./H;




A_reshape = reshape(A',[numel(s) 1]);
H_reshape = reshape(H',[numel(s) 1]);

error_freq = A_reshape .* s_fft;
error_time = ifft(error_freq);


H_norm = H ./max(abs(H),[],1) .* max(abs(s_fft_seg),[],1);
H_norm_reshape = reshape(H_norm',[numel(s) 1]);
% r_auto_correlation_norm = r_auto_correlation ./ max(r_auto_correlation, [], 1);
% instr_H(:,nn) = (instr_H(:,nn)/max(abs(instr_H(:,nn))))*max(abs(instr_st_signal_w(:,nn)));

% H = fftshift(freqz(1, a_exp));

% %%%% transfer frequency
% 
% audioOut =zeros(size(s));
% for ss = 1:num_segment
%     audioOut(ss,:) = filter([0 -a_exp(ss,2:end)], 1, s(ss,:));
% end
% 
% audioOut_reshape = reshape(audioOut',[num_segment 1]);
% 
