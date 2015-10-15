function faltning = falta(sig1, sig2)

length1 = length(sig1);
length2 = length(sig2);

nr_fft_points = length1 + length2 - 1;

spec1 = fft(sig1, nr_fft_points);
spec2 = fft(sig2, nr_fft_points);

if size(spec1) ~= size(spec2)
    spec2 = spec2.';
end

faltning = ifft(spec1.*spec2);

end