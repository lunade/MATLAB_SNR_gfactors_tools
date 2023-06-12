function  [snr_rss, snr_cov, noisecov, noisecor] = generate_cov_and_snr(meas1, meas2)
% Copyright © 2006-2012 Jonathan R. Polimeni and
%   The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
% SNR calculations reference:
%
%  Kellman P, McVeigh ER.
%  Image reconstruction in SNR units: a general method for SNR measurement.
%  Magn Reson Med. 2005 Dec;54(6):1439-47.
%  PMID: 1626157  

% have a look at the fields of the struct
meas1

% utility for viewing the data array dimensions
mrir_ice_dimensions(meas1.data)


% '1' at end is for effective noise bandwidth correction (i.e., the "0.8")
noisecov = mrir_array_stats_matrix(meas2.data, 'cov', 1);
noisecor = mrir_array_stats_matrix(meas2.data, 'cof', 1);

% this crops out in the readout direction
img = mrir_conventional_2d(meas1.data);

% assume that the image is a good approximation to the coil sensitivity profile, otherwise, smooth, fit polynomials, etc.
sens = img;

% SNR of a root-sum-of-squares combinations
snr_rss = mrir_array_SNR_rss(sens, noisecov);

% SNR of a noise covariance-weighted combination
snr_cov = mrir_array_SNR_cov(sens, noisecov);

% number of channels in the data
Ncha = mrir_ice_dimensions(meas1.data, 'cha');

% constantinides magnitude bias correction for coil array data
snr_rss_unbiased = mrir_array_combine_noise_correction(snr_rss, Ncha, 'rss');
snr_cov_unbiased = mrir_array_combine_noise_correction(snr_cov, Ncha, 'cov');


end
