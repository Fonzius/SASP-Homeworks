function[a_exp, error_time_direct,error_freq,  M, num_segment,s_fft, output_signal_time_direct, start_index, end_index ] = LPCFilter_new (audioFile, windowlength, windowtype, p )

% Import the files
[signal, fs] = audioread(audioFile);

%%%% PART 1: add windowing to the signal
%%%% segment == window 
%%%% 5 ms for SPEECH! is taken from lesson as example segment length
%%%% for harmonic and noise signal, probably we could have longer windows (
%%%% more segments in each window)
%%%% to get the window length, it need to satisfy the staionary assumption

%%%% Staionary Assumption, recording 20230317, 31:57 --> ge with time. 
%%%% we need to compute the mean and variance to make sure they don't
%%%% change in each segment.
M = windowlength;
% M = floor(5e-3*fs); % How many samples in each segment

%%%% Part 2: choose p
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


%% windowing

%%%% rectangular window
% num_segment = ceil(length(signal)/M);
% num_pad = num_segment* M -length(signal);
% paddedSignal = padarray(signal,[num_pad 0],0,'post');
% s = reshape(paddedSignal,M,num_segment)';
% s_fft = fft(paddedSignal);
% s_fft_seg = fft(s')';

%%%% hamming window -mine -wrong

win = hann(windowlength);
hop_size = floor(windowlength/4);
num_segment = ceil((length(signal)-windowlength)/hop_size + 1);
num_padding = (num_segment - 1) * hop_size + windowlength - length(signal);
signal(end+1:end+num_padding) = 0;
start_index =zeros(num_segment,1);
end_index =zeros(num_segment,1);
signal_reshape =zeros(M,num_segment);

output_signal_time_direct = zeros(length(signal),1);

for i=1:num_segment
    start_index(i) = (i-1)*hop_size + 1;
    end_index(i)= start_index(i) + windowlength - 1; 
    signal_reshape(:,i) = signal(start_index(i):end_index(i));
end
    s = win .* signal_reshape;


%   [s1,num_segment1] = windowing(signal,windowtype,windowlength);
s_fft = fft(s)';

%% auto-correlation matrix 
%%%% find a to minimize the short-time mean-squared error

r_auto_correlation = zeros(p+1, num_segment);
for ss = 1:num_segment
    for ii_kk = 0:p
        for mm = 1:M-(ii_kk)
            r_auto_correlation(ii_kk+1,ss) = r_auto_correlation(ii_kk+1,ss)+ s(mm,ss)*s(mm+ii_kk,ss);
        end
    end 
end
% wrong previous version
% M --> p+1
% auto_para_cal = zeros(M,M);
% r_auto_correlation = zeros(M);
% for ss = 1:num_segment
%     sample_segment_current = s(ss,:)';
%     % get calculation parameter matrix of auto-correlation
%     for mm = 1:M        
%         auto_para_cal(mm, 1:M-mm+1) = sample_segment_current(mm:end)'; % fill the matrix 
%     end
%     auto_correlation_current = auto_para_cal * sample_segment_current;
%     r_auto_correlation(:,ss) = auto_correlation_current';
% end

% test for normalization
% r_auto_correlation_norm = r_auto_correlation ./ max(r_auto_correlation, [], 1);


r = r_auto_correlation(2:end,:);
R = zeros(p,p,num_segment);
for ss= 1:num_segment
    for ii = 1:p
        for jj = 1:p
            R(ii,jj,ss) = r_auto_correlation(abs(ii-jj)+1,ss);
        end 
    end 
end
a = zeros(p,num_segment);
a_test = zeros(p,num_segment);
R_inverse = zeros(p,p,num_segment);

%%%%%%%%%%% VERY SLOW
for ss= 1:num_segment
    R_current_seg = R(:,:,ss);
%     R_current_seg_inv = inv(R_current_seg);
    a_test(:,ss) = R_current_seg\r(:,ss);
    R_inverse(:,:,ss) = inv(R(:,:,ss));
    a(:,ss) = R_inverse(:,:,ss) * r(:,ss);
end

r = r';
a = a';
a_exp1 = ones(size(a,1),1);
a_exp =[a_exp1 -1.*a];


% %%%% compare with lpc by matlab
aaaa = zeros(size(a_exp'));
gggg = zeros(size(a_exp'));
for ss =1:num_segment
[aaaa(:,ss),gggg(:,ss)]=lpc(s(:,ss),p);
end
aaaa =aaaa';


%% filter method by jerry

s =s';
error_time_direct = zeros(size(s));
error_freq = zeros(size(s));
error_time_test = zeros(size(s));
% synth_time_direct = zeros(size(s));
for ss =1:num_segment
   error_time_direct(ss,:) = filter(a_exp(ss,:),1,s(ss,:));
   error_freq(ss,:) =fft(error_time_direct(ss,:)')';
   error_time_test(ss,:) =ifft(error_freq(ss,:)')';
end




% %% get filter
% 
% H = zeros(size(s)); %shaping filter
% % A = zeros(size(s));
% % H_1 = zeros(size(s));
% for ss = 1:num_segment
%     [H_index,~] = freqz(1, a_exp(ss,:),"whole",M);
%     H(:,ss) = H_index';
% end
% 
% A = 1./H; %whitening filter
% 
% A_norm = normalize(A);
% H_norm = normalize(H);
% s_fft_norm = normalize(s_fft);
% % A_reshape = reshape(A,[numel(s) 1]);
% % H_reshape = reshape(H,[numel(s) 1]);
% 
% error_freq = A_norm .* s_fft_norm;
% % error_time = ifft(error_freq);
% 
% 
% % H_norm = H ./max(abs(H),[],1) .* max(abs(s_fft_seg),[],1);
% % H_norm_reshape = reshape(H_norm',[numel(s) 1]);