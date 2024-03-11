function snr = calcSNR(ref, obs)

snr = 20*log10(norm(ref)/norm(ref - obs));

end