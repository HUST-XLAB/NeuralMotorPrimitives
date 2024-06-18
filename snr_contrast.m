function snr = snr_contrast(data1, data2)
    Psingal = sum(abs(data1 .^ 2)) ./ length(data1);
    Pnoise = sum(abs(data2 .^ 2)) ./ length(data2);
    snr = 10 * log10(Psingal / Pnoise);
