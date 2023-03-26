function[H, error] = LPCFilter (audioFile)

% Import the files
[signal, fs] = audioread(audioFile);

% 5 ms is taken from lesson as example segment length
M = floor(5e-3*fs);

index = 1:M:length(signal); %resize the piano signal
s = zeros(length(index), M); %piano signal split in segments

% Create a new signal by appending zeros at the end s.t. length(paddedSignal) is
% multiple of M
paddedSignal = zeros(numel(s),1);
paddedSignal(1:length(signal)) = signal(:);

% split the signal in segments of length M
for ii = 1:length(s)
    s(ii,:) = paddedSignal(index(ii) : index(ii)+M-1);
end
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

signalAutocorr = zeros(size(s,1), size(s,2)+1); % n , p
for ii = 1:length(signalAutocorr)
    signalAutocorr(ii,:) = autocorrelation(s(ii,:),M);
end

r = signalAutocorr(:,2:end);
r1 = signalAutocorr(:,1:end-1);
% Create the symmetrical autocorrelation matrix by putting in the i,j 
% entry the i-j entry of the autocorrelation vector
R = zeros(M,M,length(s));
for ii = 1:M
    for jj = 1:M
        R(ii,jj,:) = r1(:,abs(ii-jj)+1);
    end
end

% Compute the a coefficients as R^-1 * r
a = zeros(length(s), size(s,2));
for ii = 1:length(s)
    a(ii,:) = inv(R(:,:,ii)) * r(ii,:)';
end

% Create a predicted version of the file by convolving the filter with the signal

sPredict = zeros(size(s));
for nn = 1:size(s,1)
    for kk = 1:size(s,2)
        sPredict(nn,kk) = sum(a(nn,1:kk) .* s(nn, kk:-1:1));
    end
end

predictedSignal = zeros(size(paddedSignal));
for ii = 1:size(index,2)-1
    predictedSignal(index(ii):index(ii+1)-1) = sPredict(ii,:);
end

%Z-Transoform of the predicted signal
sPredictZ = fft(predictedSignal);
N = length(sPredictZ);
z = exp(2*pi*1i/N);
k = 0:N-1;
sPredictZ = ((1/N) * sPredictZ' .* z.^(-k))';

%Z-Transform of the original signal
sZ = fft(paddedSignal);
sZ = ((1/N) * sZ' .* z.^(-k))';

%Find Z-Transform of the error
error = sPredictZ - sZ;

% Find filter H
H = (sPredictZ ./error);
