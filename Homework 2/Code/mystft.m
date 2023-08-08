function stft = mystft(signal,hop_size,window_length,padded_length)

window = hann(window_length)';
reshape_h = ceil(length(signal)/hop_size);
reshape_size = prod([window_length,reshape_h]);


% Pad with zeros to do the reshape
signal_zeros = [signal' zeros(1, abs(size(signal,1)-reshape_size))]';

% Split the original signal into overlapping segments
signal_reshape = zeros(reshape_h, window_length);

for ii = 1:reshape_h
    index = (ii-1)*hop_size+1;
    signal_reshape(ii,:) = signal_zeros(index:index+window_length-1);
end

% Apply the window to the segments
signal_window = signal_reshape.*window;

% Pad the segments to length padded_length
signal_padded = [signal_window  zeros(size(signal_window,1), padded_length-size(signal_window,2))];

% Take fft of all the segments to obtain the stft
signal_stft = fft(signal_padded, padded_length, 2);


stft = signal_stft(:, 1:padded_length/2);

end