function [fx_cont_coeffs f_dp fc_l fc_r] = fcont_gram_blend(fx, d, C, A, Q, AQ, FAQF, str)

% fftw('wisdom', str);
% load(['FC_data_d', num2str(d), '.mat']);
fourPts = length(fx) + C;

fr = fx((length(fx) - (d-1)):length(fx));
fl = fx(1:d);
% f_dp = zeros(fourPts, 1);

% fc_r = A*(Q.'*fr);
% fc_l = flipud(A*(Q.'*flipud(fl)));
fc_r = AQ*fr;
fc_l = FAQF*fl;

% f_dp(1 : length(fx)) = fx;
% f_dp(length(fx) + 1 : end) = fc_l + fc_r;
f_dp = [fx; fc_l + fc_r];

fx_cont_coeffs = fft(f_dp)/fourPts;
% fx_cont_coeffs = fft(f_dp);
return
