function [signal] = removeZeros(signal)

for ii = size(signal):-1:find(signal,1,'last')
    signal(ii) = [];
end

for ii = 1:find(signal,1)-1
    signal(1) = [];
end

end