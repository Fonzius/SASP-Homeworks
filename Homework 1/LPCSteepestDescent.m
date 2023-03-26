function[H] = LPCSteepestDescent(audioFile)

audioFile = "piano.wav";
[signal, fs] = audioread(audioFile);

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

R = zeros(M,M,length(s));
for ii = 1:length(s)
    R(:,:,ii) = corrmtx(s(ii,:),size(s,2)-1);
end

mu = 1/max(eig(R));

tau = 5 * ceil(1/(2*mu+min(eig(R)))); %CHECK HOW MANY TAU WE NEED TO CONVERGE TO SOMETHING NICE

w = zeros(length(s),M);
u = zeros(length(s),M);


for tt = 1:tau
    for ii = 1:length(s)
        temp = conv(s(ii,:),w(ii,:));
        u(ii,:) = temp(ceil(length(temp)/2):end);
    end
    
    p = zeros(size(s));
    for ii = 1:lenght(s)
        for jj = 1:M
            p(ii,jj) = avg(u(ii,1:jj) .* conj(s(ii,1:jj)));
        end
    end
    
    
    
    for ii = 1:length(s)
        w(ii,:) = w(ii,:) * mu*(p(ii,:)-R(:,:,ii)*w(ii,:)');
    end
end




