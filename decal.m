function [res] = decal(Sdec, dec, signal)

res = zeros(length(signal),1);
for i=1:length(signal)-1
    if ( i+dec >0 && i+dec <length(Sdec))
            res(i) = Sdec(i+dec);
    end
end