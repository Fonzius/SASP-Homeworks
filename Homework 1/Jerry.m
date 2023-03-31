clc; clearvars; close all;

[speech, ~] = audioread("./speech.wav");
[music, Fs] = audioread("./piano.wav");

p_speech = 32;
p_music = 128;
n_steps = 2000;
fac = 0.35;

L = length(speech);

win_length = 2048;
win = hann(win_length);
hop_size = floor(win_length/4);

N_windows = ceil((L-win_length)/hop_size + 1);
num_zeroes_to_fill = (N_windows - 1) * hop_size + win_length - L;

speech(end+1:end+num_zeroes_to_fill) = 0;
music(end+1:end+num_zeroes_to_fill) = 0;
L = length(speech);

output_close_form = zeros(L, 1);
output_steepest = zeros(L, 1);


for i=3:N_windows
    
    start_index = (i-1)*hop_size + 1;
    end_index = start_index + win_length - 1;

    music_segment = win .* music(start_index:end_index, 1);
    speech_segment = win .* speech(start_index:end_index, 1);
    
    [r_music, lags_music] = xcorr(music_segment);
    [r_speech, lags_speech] = xcorr(speech_segment);
    
    r_music = r_music(lags_music>=0);
    r_speech = r_speech(lags_speech>=0);

    R_music = toeplitz(r_music(1:p_music));
    R_speech = toeplitz(r_speech(1:p_speech));

    r_music = r_music(2:p_music+1);
    r_speech = r_speech(2:p_speech+1);
    
    % close form
    a_music = R_music\r_music;
    a_speech = R_speech\r_speech;

    % steepest
    eigenvalues_speech = eig(R_speech);
    eigenvalues_music = eig(R_music);
    mu_speech = fac * 2/max(eigenvalues_speech);
    mu_music = fac * 2/max(eigenvalues_music);
    
    w_speech = zeros(p_speech,1);
    w_music = zeros(p_music,1);
    for n = 1:n_steps
        w_speech = w_speech + mu_speech*(r_speech-R_speech*w_speech);
        w_music = w_music + mu_music*(r_music-R_music*w_music);
    end

    A_music_close_form = [1; -a_music];
    A_speech_close_form = [1; -a_speech];
    A_music_steepest = [1; -w_music];
    A_speech_steepest = [1; -w_speech];

    music_error_close_form = filter(A_music_close_form, 1, music_segment);
    speech_error_close_form = filter(A_speech_close_form, 1, speech_segment);
    music_error_steepest = filter(A_music_steepest, 1, music_segment);
    audio_synth_close_form = filter(1, A_speech_close_form, music_error_close_form);
    audio_synth_steepest = filter(1, A_speech_steepest, music_error_steepest);

    %%% go to frequency domain
    speech_segment_freq = fft(speech_segment);
    speech_errorfreq = fft(speech_error_close_form);
    music_errorfreq = fft(music_error_close_form);
    Hspeech_freq = speech_segment_freq./ speech_errorfreq ;
    talking_freq = music_errorfreq.* Hspeech_freq;
    talking_time = ifft(talking_freq);
    %%%% check if audio_synth_close_form ==  talking_time
   
  
    output_close_form(start_index:end_index, 1) = output_close_form(start_index:end_index, 1) + audio_synth_close_form;
    output_steepest(start_index:end_index, 1) = output_steepest(start_index:end_index, 1) + audio_synth_steepest;
end

output_close_form = output_close_form / max(abs(output_close_form));
output_steepest = output_steepest / max(abs(output_steepest));

audiowrite('output_close_form.wav', output_close_form, Fs)
sound(output_close_form,44100)
audiowrite('output_steepest.wav', output_steepest, Fs)

