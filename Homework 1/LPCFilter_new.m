function[a_exp, fs, M, num_segment,s_fft, ...
            length_tot_signal, start_index, end_index ] = ...
                LPCFilter_new (audioFile, windowlength, p , method )

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

win = hann(windowlength);
hop_size = floor(windowlength/4);
num_segment = ceil((length(signal)-windowlength)/hop_size + 1);
num_padding = (num_segment - 1) * hop_size + windowlength - length(signal);
signal(end+1:end+num_padding) = 0;
start_index =zeros(num_segment,1);
end_index =zeros(num_segment,1);
signal_reshape =zeros(M,num_segment);

length_tot_signal = length(signal);


for i=1:num_segment
    start_index(i) = (i-1)*hop_size + 1;
    end_index(i)= start_index(i) + windowlength - 1; 
    signal_reshape(:,i) = signal(start_index(i):end_index(i));
end
    s = win .* signal_reshape;
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
%% Close form method
if method == 1

    a = zeros(p,num_segment);
%     a_test = zeros(p,num_segment);
    R_inverse = zeros(p,p,num_segment);
    
    %%%%%%%%%%% VERY SLOW
    for ss= 1:num_segment
%         R_current_seg = R(:,:,ss);
    %     R_current_seg_inv = inv(R_current_seg);
%         a_test(:,ss) = R_current_seg\r(:,ss);
        R_inverse(:,:,ss) = inv(R(:,:,ss));
        a(:,ss) = R_inverse(:,:,ss) * r(:,ss);
    end    
%     a = a';
%     a_exp1 = ones(size(a,1),1);
%     a_exp =[a_exp1 -1.*a];    
%     % %%%% compare with lpc by matlab
%     aaaa = zeros(size(a_exp'));
%     gggg = zeros(size(a_exp'));
%     for ss =1:num_segment
%     [aaaa(:,ss),gggg(:,ss)]=lpc(s(:,ss),p);
%     end
%     aaaa =aaaa';
%     gggg =gggg';


%% Steepest decent method
elseif method == 2
    n_steps = 2000;

    lambda = zeros(p,num_segment); %check dimenstion!
    
    a = zeros(p,num_segment);

    for ss= 1:num_segment
        lambda(:,ss) = eig(R(:,:,ss));
    end
    % The necessary and sufficient condition for the convergence or stability
    % 0 < lambda <2/lambdaMax
    lambdaMax = max(lambda);
    lambdaMin = min(lambda);

    fac = 0.5;
    mu = fac .* 2 ./ lambdaMax;
    tau = 1./(2.*mu.*lambdaMin);
%     a_loop = zeros(p,num_segment,n_steps+1);

    for ss= 1:num_segment
        for n = 1:n_steps            
            a(:,ss) = a(:,ss) + mu(1,ss).*(r(:,ss)-R(:,:,ss)*a(:,ss)); 
%             a_loop(:,ss,n+1) = a(:,ss);          
        end    
    end


    if test_transient ==1
        a = zeros(p,num_segment);
        a_loop = zeros(p,num_segment,n_steps+1);
        for ss= 1:num_segment
        for n = 1:n_steps            
            a(:,ss) = a(:,ss) + mu(1,ss).*(r(:,ss)-R(:,:,ss)*a(:,ss)); 
            a_loop(:,ss,n+1) = a(:,ss);          
        end    
    end
%     a_loop_exp = -1.* a_loop;
%     a_new_row = ones(1, num_segment, n_steps+1); 
%     a_loop_exp = cat(1, a_new_row, a_loop_exp); % Concatenate the two arrays along the first dimension

    %%%%% for Transient behaviour test

    a_wen = zeros(p,num_segment);
    R_inverse = zeros(p,p,num_segment);
    %%%%%%%%%%% VERY SLOW
    for ss= 1:num_segment
        R_inverse(:,:,ss) = inv(R(:,:,ss));
        a_wen(:,ss) = R_inverse(:,:,ss) * r(:,ss);
    end 

%     a_lpc = zeros(p+1,num_segment);
%     for ss =1:num_segment
%     [acurrent,~] = lpc(s(:,ss),p);   
%      a_lpc(:,ss) = acurrent';
%     end
    % weight-error vector
    c = a_loop - repmat(a_wen, [1, 1, size(a_loop, 3)]);
    Q = zeros(p,p,num_segment);
    QH = zeros(p,p,num_segment);
    for ss= 1:num_segment
        [Q(:,:,ss), ~] = eig(R(:,:,ss));
        QH(:,:,ss) =conj(Q(:,:,ss)');
    end
    v =zeros(size(c));
    for ss= 1:num_segment
        v(:,:,ss) = QH(:,:,ss) * c(:,:,ss);
    end

    %%v: 10(eigenvalue)*1656(seg)*2001(stepinerate)
    figure()
    vcur = reshape(v(1,1,:), [1, 2001]);
    plot(1:2001, vcur)
    end
    
else 
    error('Invalid method!!!')
end

a = a';
a_exp1 = ones(size(a,1),1);
a_exp =[a_exp1 -1.*a];   

