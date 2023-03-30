clc; close all; clear all;

%==========================================================================
%                           DAAP HW1
%                           
%==========================================================================

% xcorr, circshift, cell, window

instr_name = "piano.wav";
instr_name_recon = "piano_recon.wav";
speech_name = "speech.wav";
st_res_name = "whitening_filter_res.wav";
st_res_name_test = "whitening_filter_res_test.wav";
st_res_sub_name = "subtracted_res.wav";
st_res_sub_name_test = "subtracted_res.wav";


[instr_t,instr_Fs] = audioread(instr_name);
[speech_t,speech_Fs] = audioread(speech_name);

instr_t_len =length(instr_t);

taps = 700;
p = taps;
wl =  1024;  % windowLenght =

%t_buckets = ceil(instr_t_len/wl);
%instr_t(end+1:wl*t_buckets) = 0;
%st_signal = reshape(instr_t,[wl t_buckets]);
[st_signal,chunksNum] = windowing(instr_t,"hamming",1024);


t_buckets = chunksNum;
r = zeros(p+1, t_buckets);

M = wl;
for nn = 1:t_buckets
    for ii_kk = 0:p
        for mm = 1:M-(ii_kk)
            r(ii_kk+1,nn) = r(ii_kk+1,nn)+ st_signal(mm,nn)*st_signal(mm+ii_kk,nn);
        end
    end 
end
r_0 = r(1:end-1,:);
figure
plot(r_0);
r_1 = r(2:end,:);
R = zeros(p,p,t_buckets);
a = zeros(p, t_buckets);
for nn = 1:t_buckets
    R(:,:,nn) = toeplitz(r_0(:,nn));
    a(:,nn) = R(:,:,nn)\r_1(:,nn);
end


R_1 = R(:,:,1);
%% Reconstruction of the prediction using prediction filters 

%reconstructed_Sig = zeros(wl, t_buckets);

%for nn = 1:t_buckets
%    for mm = p+1:M
%        for kk = 1:p
%            reconstructed_Sig(mm,nn) = reconstructed_Sig(mm,nn) ...
%                + a(kk,nn)*st_signal(mm-kk+1,nn);
%        end
%    end
%end
%reconstructed_Sig = reshape(reconstructed_Sig,[t_buckets *wl 1]);
%reconstructed_Sig = adding(st_)
%audiowrite(instr_name_recon,reconstructed_Sig,instr_Fs);

%% Z_transform of the filter represented as spectrum of the filter
%n = wl;

%w = linspace(-instr_Fs/2,instr_Fs/2,n);

%A = zeros(length(w),t_buckets);
%z= exp(1i*w);
%sum_A =zeros(length(w),t_buckets);
%for nn=1:t_buckets
%    for kk = 1:p
%        sum_A(:,nn) = a(kk,nn)*z.^(-kk);
%    end
%end
%A = 1- sum_A;
%Aishift = zeros(wl,t_buckets);
%for nn = 1: t_buckets
%    Aishift(:,nn) = ifftshift(A(:,nn));
%end 
%A1 = Aishift(1:512,:);

%% STFT of the signal 
w1   = linspace(0,instr_Fs/2,wl/2);
st_signal_w = zeros(wl,t_buckets);
for nn = 1:t_buckets
    st_signal_w(:,nn) = fft(st_signal(:,nn));
end

Fs = instr_Fs;
figure
t =linspace(0,instr_t_len/instr_Fs,t_buckets);
st_signal_w1 = st_signal_w(1:512,:);
surf(t,w1,abs(st_signal_w1),EdgeColor="none");
ylim([0 5000])
view(0,90)
disp(st_signal_w(1,400))


%% Z tranform other method
a_ones = ones(1,t_buckets);
a_1 = zeros(p+1,t_buckets);
A = zeros(wl,t_buckets);
H = zeros(wl,t_buckets);

for nn = 1: t_buckets
    a_1(:,nn) = vertcat(1, -a(:,nn));
    
    [A(:,nn),~] = freqz(a_1(:,nn),1,"whole",wl);
    A(:,nn) = A(:,nn)/(mean(A(:,nn)/mean(st_signal_w(:,nn))));
    [H(:,nn),w] = freqz(1,a_1(:,nn),"whole",wl);
    H(:,nn) = H(:,nn)/(mean(H(:,nn))/mean(st_signal_w(:,nn)));
    %A(:,nn) = 1/(H(:,nn));
    
end 
figure
plot(w/pi,A(:,800));
hold on;
plot(w/pi,H(:,800));


%% piano signal filtering
st_signal_w = zeros(wl,t_buckets);
for nn = 1:t_buckets
    st_signal_w(:,nn) = fft(st_signal(:,nn));
end


st_res =  zeros(wl, t_buckets); 
st_res_w = zeros(wl,t_buckets);
for nn = 1:t_buckets
    st_res_w(:,nn) =  st_signal_w(:,nn).* A(:,nn);
    %st_res_w(:,nn) = st_signal_w(:,nn)./ H(:,nn);
    st_res(:,nn) = ifft(st_res_w(:,nn));
end

%st_res_lin = reshape(st_res,[t_buckets *wl 1]);
st_res_lin = adding(st_res,0.5,1024);

audiowrite(st_res_name,abs(st_res_lin),instr_Fs);

instr_t(end+1:length(st_res_lin)) = 0;
st_res_form = instr_t- st_res_lin;

audiowrite(st_res_sub_name,abs(st_res_form),instr_Fs);





%% test using matlab lpc

H_test = zeros(wl,t_buckets);
for nn = 1:t_buckets
   a_test(:,nn) = lpc(st_signal(:,nn), p);
   %a_1_test(:,nn) = vertcat(1, -a_test(:,nn));
    
    [A_test(:,nn),~] = freqz(a_test(:,nn),1,"whole",wl);
    %A_test(:,nn) = A_test(:,nn)/(mean(A_test(:,nn)/mean(st_signal_w(:,nn))));
    [H_test(:,nn),w] = freqz(1,a_test(:,nn),"whole",wl);
    %H_test(:,nn) = H_test(:,nn)/(mean(H_test(:,nn))/mean(st_signal_w(:,nn)));
    %A(:,nn) = 1/(H(:,nn));
end

%% plot of a random element of the matrix for the signal and for the filter 
figure
plot(w,10*log10(abs(st_signal_w(:,600)).^2))
%plot(w,abs(st_signal_w(:,600)))

hold on 
plot(w,10*log10(abs(H(:,600)).^2))
%plot(w,abs(H(:,600)))
hold on 
plot(w,10*log10(abs(H_test(:,600)).^2),"LineStyle","--")
%plot(w,abs(H(:,600)))

xlim([0,1])
legend("signal chunk","filter shape","test filter shape")

%% piano signal filtering test


st_res_test =  zeros(wl, t_buckets); 
st_res_w_test = zeros(wl,nn);
for nn = 1:t_buckets
    st_res_w_test(:,nn) =  st_signal_w(:,nn).* A_test(:,nn);
    %st_res_w(:,nn) = st_signal_w(:,nn)./ H(:,nn);
    st_res_test(:,nn) = ifft(st_res_w(:,nn));
end

%st_res_lin = reshape(st_res,[t_buckets *wl 1]);
st_res_lin_test = adding(st_res,0.5,1024);

audiowrite(st_res_name_test,abs(st_res_lin),instr_Fs);

instr_t(end+1:length(st_res_lin_test)) = 0;
st_res_form_test = instr_t- st_res_lin_test;

audiowrite(st_res_sub_name_test,abs(st_res_form_test),instr_Fs);