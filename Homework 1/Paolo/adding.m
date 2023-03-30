function [addedSignal] = adding(shortTimeSignal,windowOverLap,windowLength)
% =========================================================================
% windowOverLap: overLap of the windows as a fractional number
% 
%
% =========================================================================
R = windowOverLap * windowLength;
M = windowLength;
[windowLength,chunks_num] = size(shortTimeSignal);
signalLength = (chunks_num * R)+(M-R);
addedSignal = zeros(signalLength,1);
ii=1;
for nn = 1:R:(signalLength-(M-R))
    s_chunk = shortTimeSignal(:,ii);
    addedSignal( nn:nn+(M-1),1) = addedSignal(nn:nn+(M-1),1) + s_chunk;
    ii =ii+1;
    %disp(strcat("nn=", num2str(nn)))
end

%audiowrite("addedSignal.wav",addedSignal,44100);