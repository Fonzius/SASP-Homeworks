function [a] = autocorrelation(signal, max_lag)

% Compute the mean value of x and subtract it to obtain a zero-mean signal
mean_signal = mean(signal);
signal_zero_mean = signal - mean_signal;

% Create a vector of lags and initialize the autocorrelation vector
lags = 0:max_lag;
a = zeros(1,max_lag+1);

% Compute the autocorrelation for each lag
for k = 1:length(lags)
    lag = lags(k);
    a(k) = sum(signal_zero_mean(1:end-lag) .* signal_zero_mean(1+lag:end));
end

% Normalize the autocorrelation values
N = length(signal);
var_signal = var(signal);
a = a / (N * var_signal);