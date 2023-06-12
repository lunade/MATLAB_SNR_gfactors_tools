% function dcm2niix_wrapper(dicom_folder, nifti_folder, nifti_filename)
%
% nifti filename does not need extension .nii


function dcm2niix_wrapper(dicom_folder, nifti_folder, nifti_filename)
%
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

nifti_folder
nifti_filename
dicom_folder


unix(['/usr/pubsw/packages/mricrogl/new/dcm2niix -o ',[dicom_folder,'/' nifti_folder],' -f ', nifti_filename ' ', dicom_folder]);

% ['dcm2niix -w 1 -o ', nifti_folder, ' -f ',nifti_filename, ' ', dicom_folder]



% unix(['/usr/pubsw/packages/mricrogl/new/dcm2niix -w 1 -o ', nifti_folder, ' -f ',nifti_filename, ' ', dicom_folder]);

 
 
 
 
