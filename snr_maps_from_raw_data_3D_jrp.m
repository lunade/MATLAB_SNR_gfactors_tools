function [snr_cov, snr_rss, img, noisecov, noisecorr, kspace] = snr_maps_from_raw_data_3D_jrp(sig_path, noise_path,fix_slices_switch);
% function [snr_cov, snr_rss, img, noisecov, noisecorr, kspace] = snr_maps_from_raw_data_3D_jrp(sig_path, noise_path,fix_slices_switch); 
%Copyright © 2006-2012 Jonathan R. Polimeni and
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
% jonathan polimeni <jonp@nmr.mgh.harvard.edu>
% Friday, April 29, 2011 15:22:38 -0400 



which read_meas_dat
% read_meas_dat is compatible with VB* baselines, NOT VD*!
meas1 = read_meas_dat(sig_path);


if isfield(meas1,'noiseadjscan') == 0

    meas2 = read_meas_dat(noise_path);
else
    meas2.data = meas1.noiseadjscan;
end

% have a look at the fields of the struct
meas1

% utility for viewing the data array dimensions
mrir_ice_dimensions(meas1.data)

% downsample the readout direction
meas1.data=meas1.data(1:2:end,:,:,:,:,:,:,:,:,:,:);
meas2.data=meas2.data(1:2:end,:,:,:,:,:,:,:,:,:,:);

% '1' at end is for effective noise bandwidth correction (i.e., the "0.8")
noisecov = mrir_array_stats_matrix(meas2.data, 'cov', 1);

noisecorr = mrir_array_stats_matrix(meas2.data, 'cof', 1);



% this crops out in the readout direction
% img = mrir_conventional_3d(meas1.data);
temp=permute(meas1.data,[1 2 9 3 4 5 6 7 8]);  
for cc=1:numel(temp(1,1,1,:,1,1,1,1,1,1))
    temp2(:,:,:,cc) = fftshift(fftn(fftshift(temp(:,:,:,cc,1,1,1,1,1,1))));
end

% figure(2020),imagesc(vol2mos(abs(temp2(:,:,:,20))))
img=permute(temp2,[1 2 4 5 6 7 8 9 3]);


% assume that the image is a good approximation to the coil sensitivity profile, otherwise, smooth, fit polynomials, etc.
sens = img;

% SNR of a root-sum-of-squares combinations
snr_rss = squeeze(mrir_array_SNR_rss(sens, noisecov));

% SNR of a noise covariance-weighted combination
snr_cov = squeeze(mrir_array_SNR_cov(sens, noisecov));

% number of channels in the data
Ncha = mrir_ice_dimensions(meas1.data, 'cha');

% constantinides magnitude bias correction for coil array data
% snr_rss_unbiased = mrir_array_combine_noise_correction(snr_rss, Ncha, 'rss');
% snr_cov_unbiased = mrir_array_combine_noise_correction(snr_cov, Ncha, 'cov');
if mod(numel(snr_cov(1,1,:)),0) == 1
    snr_cov=snr_cov(:,:,1:end-1);
    snr_rss=snr_rss(:,:,1:end-1);
end

if fix_slices_switch == 1
    snr_cov = fix_slices(snr_cov);
    snr_rss = fix_slices(snr_rss);
end

% fix slice interleaves


img=(permute(squeeze(img),[1 2 4 3]));

for ss=1:numel(img(1,1,:,1))
    for cc=1:numel(img(1,1,1,:))
        kspace(:,:,ss,cc)=fftshift(fft2(fftshift(img(:,:,ss,cc))));
    end
end





% use BART to estimate coil sensitivities for g-factor calculation
%   sensitivities = bart('ecalib', kspace);
%   image_out = bart('pics -l1 -r0.001', kspace, sensitivities);

%  calib = bart('ecalib -r 24 -m 1 -d 7 -c 0.8 ', fft3c(img_3d));
%  sens=single(calib(:,:,:,:,1));  

 
 
% for cc=1:numel(img(1,1,1,:))
%     img_temp(:,:,:,cc) = fix_slices(img(:,:,:,cc));
% end
% img=img_temp;

end

 function out = fix_slices(in)
    if mod(numel(in(1,1,:)),2) == 0
        out=zeros(size(in));
        out(:,:,1:2:end) = in(:,:,end/2+1:end);
        out(:,:,2:2:end) = in(:,:,1:end/2);
    else  % odd number of slices
        out=zeros(size(in));
        
        out(:,:,1:2:end) = in(:,:,1:floor(end/2)+1);
        out(:,:,2:2:end) = in(:,:,1:floor(end/2));
    end
        out=rot90(out);

 end



% say we have image data AND coil sensitivity profile data:
% sens = some_function(some_data);
% 
% % sens and img are COMPLEX-VALUED
% snr_opt = mrir_array_SNR_opt(sens, noisecov, img);
% 
% snr_opt_unbiased = mrir_array_combine_noise_correction(snr_cov, Ncha, 'opt');
% 
% % this is equivalent to:
% snr_opt_unbiased = mrir_array_combine_noise_correction(snr_cov, 1, 'opt');


