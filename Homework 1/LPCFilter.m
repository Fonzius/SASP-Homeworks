<<<<<<< Updated upstream
function[s_predict] = LPCFilter (audioFile)
=======
function[a] = LPCFilter (audioFile)
>>>>>>> Stashed changes

% Import the files
[signal, fs] = audioread(audioFile);

% 5 ms is taken from lesson as example segment length
M = floor(5e-3*fs); % How many samples in each segment


%%%%//// Method by Xinmeng
num_segment = ceil(length(signal)/M);
num_pad = num_segment* M -length(signal);
paddedSignal = padarray(signal,[num_pad 0],0,'post');
s = reshape(paddedSignal,M,num_segment)';

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
a = a'; % a vector
r_auto_correlation = r_auto_correlation';
<<<<<<< Updated upstream
% R = permute(R, [2 1 3]);


%%%%//////////////////////////////////


%%
% Create a predicted version of the file by convolving the filter with the signal
=======
>>>>>>> Stashed changes


a_exp1 = ones(size(a,1),1);
a_exp =[a_exp1 -1.*a];
H = zeros(size(a_exp)); % shaping filter H
for ii = 1:size(a_exp,1)
    H(ii,:) = freqz(1, a_exp(ii,:), M);
end
A = 1./H; % whitening filter A

<<<<<<< Updated upstream
ww=1;


%%%%/// predicted version by Marco
% sPredict = zeros(size(s));
% for nn = 1:size(s,1)
%     for kk = 1:size(s,2)
%         sPredict(nn,kk) = sum(a(nn,1:kk) .* s(nn, kk:-1:1));
%     end
% end
% 
% predictedSignal = zeros(size(paddedSignal));
% for ii = 1:size(index,2)-1
%     predictedSignal(index(ii):index(ii+1)-1) = sPredict(ii,:);
% end


%%%%%%%%%%%%%%%
% % Z-Transform of the predicted signal by Xinmeng
% 
% sPredict_fft = zeros(size(s_predict)); % 用于存储每行 FFT 的结果
% 
% for i = 1:num_segment
%     sPredict_fft(i, :) = fft(s_predict(i, :));
% end
% 
% % Z-Transoform of the predicted signal by Marco
% sPredictZ = fft(predictedSignal);
% N = length(sPredictZ);
% z = exp(2*pi*1i/N);
% k = 0:N-1;
% sPredictZ = ((1/N) * sPredictZ' .* z.^(-k))';
% 
% %Z-Transform of the original signal
% sZ = fft(paddedSignal);
% sZ = ((1/N) * sZ' .* z.^(-k))';
% 
% %Find Z-Transform of the error
% error = sPredictZ - sZ;
% 
% % Find filter H
% H = (sPredictZ ./error);
=======
a_exp0 = zeros(size(a,1),1);
a_exp_p =[a_exp0 a];

P = zeros(size(a_exp_p));
for ii = 1:size(a_exp_p,1)
    P(ii,:) = freqz(a_exp_p(ii,:),1,M);
end


piano_fft = fft(s');
piano_fft = piano_fft';

error = A.*piano_fft;
error_reshape = reshape(error',[848980 1]);
error_time = ifft(error_reshape')';

piano_hat = P .* piano_fft;
piano_hat_reshape = reshape(piano_hat',[848980 1]);




>>>>>>> Stashed changes
