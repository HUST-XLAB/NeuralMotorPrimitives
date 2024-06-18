function snr_coef = calculate_series_snr(data_series, windowWidth)
    n = length(data_series);
    Psingal = [];

    for i = 1:(n - windowWidth + 1)
        Pi = sum(abs(data_series(i:i + windowWidth - 1)) .^ 2) ./ windowWidth;
        Psingal = cat(2, Psingal, Pi);
    end

    snr_coef = 10 * log10(max(Psingal) / min(Psingal(350:400)));
