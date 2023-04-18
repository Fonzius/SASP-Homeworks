% HW1 - DAAP
% by Xinmeng Luan, Marco Bernasconi
% 31 Mar,2023

function[a_exp, fs, M, num_segment,s_fft, ...
    length_tot_signal, start_index, end_index ] = ...
                LPCFilter (audioFile, windowlength, p , method, test_transient )
if method == 1
    test_transient = 0;
end

[signal, fs] = audioread(audioFile);
M = windowlength;


%% Windowing

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
    % find a to minimize the short-time mean-squared error   
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
    
    % inverse matrix is slow
    for ss= 1:num_segment
        R_inverse(:,:,ss) = inv(R(:,:,ss));
        a(:,ss) = R_inverse(:,:,ss) * r(:,ss);
    end    
   
    % compare with lpc by matlab
    % aaaa = zeros(size(a_exp'));
    % gggg = zeros(size(a_exp'));
    % for ss =1:num_segment
    % [aaaa(:,ss),gggg(:,ss)]=lpc(s(:,ss),p);
    % end
    % aaaa =aaaa';
    % gggg =gggg';


%% Steepest decent method
elseif method == 2
    n_steps = 20000;

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

    for ss= 1:num_segment
        for n = 1:n_steps            
            a(:,ss) = a(:,ss) + mu(1,ss).*(r(:,ss)-R(:,:,ss)*a(:,ss));        
        end    
    end

   %%%%% for Transient behaviour test
    if test_transient ==1
        a = zeros(p,num_segment);
        a_loop = zeros(p,num_segment,n_steps+1);
        for ss= 1:num_segment
            for n = 1:n_steps            
                a(:,ss) = a(:,ss) + mu(1,ss).*(r(:,ss)-R(:,:,ss)*a(:,ss)); 
                a_loop(:,ss,n+1) = a(:,ss);          
            end    
        end
       
        a_wen = zeros(p,num_segment);
        R_inverse = zeros(p,p,num_segment);

        for ss= 1:num_segment
            R_inverse(:,:,ss) = inv(R(:,:,ss));
            a_wen(:,ss) = R_inverse(:,:,ss) * r(:,ss);
        end 

        % weight-error vector
        c = a_loop - repmat(a_wen, [1, 1, size(a_loop, 3)]);
        figure()
         for ii = 1:size(c,1)
            for jj =1:size(c,2)
                plot(1:size(c,3),reshape(c(ii,jj,:), [1, size(c,3)]));
                title(['Weight-error vector    k=' num2str(ii) ',segment=' num2str(jj)])
                pause(0.0000001);
            end
         end
         Q = zeros(p,p,num_segment);

        for ss = 1:num_segment
            [V, ~] = eig(R(:, :, ss));  
            Q(:, :, ss) = V;            
            for i = 1:size(Q, 2)
                Q(:, i, ss) = Q(:, i, ss) / norm(Q(:, i, ss));  
            end
        end
    
        QH = conj(permute(Q, [2 1 3]));   
        v =zeros(size(c));
        for ss= 1:num_segment
            v(:,:,ss) = QH(:,:,ss) * c(:,:,ss);
        end    
      
        figure()
        for ii = 1:size(v,1)
            for jj =1:size(v,2)
                vcur = reshape(v(ii,jj,:), [1, size(v,3)]);
                plot(1:size(v,3), abs(vcur))
                xlim([0 20])
                title(['Transient behaviour    k=' num2str(ii) ',segment=' num2str(jj)])
                pause(0.0000001); 
            end
        end
    end
    
else 
    error('Invalid method!!!')
end

a = a';
a_exp1 = ones(size(a,1),1);
a_exp =[a_exp1 -1.*a];   

