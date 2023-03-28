function[a] = LPCFilter (audioFile)

% Import the files
[signal, fs] = audioread(audioFile);

% 5 ms is taken from lesson as example segment length
M = floor(5e-3*fs); % How many samples in each segment

%%%%//// Method by Xinmeng
num_segment = ceil(length(signal)/M);
num_pad = num_segment* M -length(signal);
paddedSignal = padarray(signal,[num_pad 0],0,'post');
s = reshape(paddedSignal,M,num_segment)';

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

% 
% r_auto_correlation = r_auto_correlation';
% % R = permute(R, [2 1 3]);
% 
% % Create a predicted version of the file by convolving the filter with the signal
% 
% %%%%/// predicted signal by Xinmeng
% s_predict = zeros(num_segment,M-1);
% s_predict_para_cal = zeros(M-1,M-1);
% 
% for ss = 1:num_segment
%     sample_segment_current = s(ss,1:end-1)';
%     % get calculation parameter matrix of predict signal
%     for mm = 1:M-1        
%         s_predict_para_cal(mm:M-1, mm) = sample_segment_current(1:M-mm)'; % fill the matrix 
%     end
%     s_predict_current = s_predict_para_cal * sample_segment_current;
%     s_predict(ss,:) = s_predict_current';
% end
% 
% A = zeros(size(s));
% for ii = 1:size(s,1)
%     A(ii,:) = freqz(1, a(ii,:), M);
% end
% 
% whitenedSignal = zeros(size(s));
% for ii = 1:size(s,1)
%     whitenedSignal(ii,:) = fft(s(ii,:));
%     whitenedSignal(ii,:) = whitenedSignal(ii,:) .* A(ii,:);
% end
% 
% %whitenedSignal = reshape(whitenedSignal, size(paddedSignal));
% 
% H = A;