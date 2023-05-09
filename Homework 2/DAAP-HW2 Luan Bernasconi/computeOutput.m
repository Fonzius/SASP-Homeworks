function computeOutput(labels,y1_stft,padded_length,y1,reshape_h,hop_size,fs)

mask1 = labels;
mask1(mask1 == 2 | mask1 == 3 ) = 0;
mask1(mask1 == 1) = 1;

mask1_mirrored = [mask1, flip(mask1,2)];

mask2 = labels;
mask2(mask2 == 1 | mask2 == 3 ) = 0;
mask2(mask2 == 2) = 1;

mask2_mirrored = [mask2, flip(mask2,2)];

mask3 = labels;
mask3(mask3 == 1 | mask3 == 2 ) = 0;
mask3(mask3 == 3) = 1;

mask3_mirrored = [mask3, flip(mask3,2)];


% Apply the filtering with the mask
f1_stft = y1_stft.*mask1_mirrored;
f2_stft = y1_stft.*mask2_mirrored;
f3_stft = y1_stft.*mask3_mirrored;

% Take the inverse fft to obtain the segments
f1_reshape = ifft(f1_stft, padded_length, 2);
f2_reshape = ifft(f2_stft, padded_length, 2);
f3_reshape = ifft(f3_stft, padded_length, 2);

% Sum all the segments back together with overlap to obtain 
% the resulting signal
f1 = zeros(length(y1)+2*padded_length, 1);
f2 = zeros(length(y1)+2*padded_length, 1);
f3 = zeros(length(y1)+2*padded_length, 1);


for ii = 1:reshape_h
    index = ii*hop_size;
    f1(index:index+padded_length-1) = f1(index:index+padded_length-1) + f1_reshape(ii,:)';
    f2(index:index+padded_length-1) = f2(index:index+padded_length-1) + f2_reshape(ii,:)';
    f3(index:index+padded_length-1) = f3(index:index+padded_length-1) + f3_reshape(ii,:)';
end


f1_true = audioread("DAAP_HW2_reference_files/s1.wav");
f2_true = audioread("DAAP_HW2_reference_files/s2.wav");
f3_true = audioread("DAAP_HW2_reference_files/s3.wav");

f = [f1,f2,f3];
f_true = [f1_true, f2_true, f3_true];

maxCorr = zeros(3,3);
for ii = 1:3
    for jj = 1:3
        maxCorr(ii,jj) = max(abs((xcorr(f(:,ii),f_true(:,jj)))));
    end
end
[M,MaxIndexes] = max(maxCorr,[],2);

f1_reconstructed = removeZeros(f(:,MaxIndexes(1)));
f2_reconstructed = removeZeros(f(:,MaxIndexes(2)));
f3_reconstructed = removeZeros(f(:,MaxIndexes(3)));


audiowrite("Luan_Bernasconi_extraFeature_f1.wav",real(f1_reconstructed), fs);
audiowrite("Luan_Bernasconi_extraFeature_f2.wav",real(f2_reconstructed), fs);
audiowrite("Luan_Bernasconi_extraFeature_f3.wav",real(f3_reconstructed), fs);
end