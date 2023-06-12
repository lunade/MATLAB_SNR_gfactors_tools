function [snr_meas, snr_true, std_ratio, avg_ratio, snr_bias] ...
    = mrir_array_combine_noise_correction_lookup(snr_true, Ncha, varargin)
%MRIR_ARRAY_COMBINE_NOISE_CORRECTION_LOOKUP
%
% [snr_meas, snr_true] = mrir_array_combine_noise_correction_lookup(snr_true, Ncha)
%
%
% example:
%
%   mrir_array_combine_noise_correction_lookup(0, 1);


% Henkelman RM.
% Measurement of signal intensities in the presence of noise in MR images.
% Med Phys. 1985 Mar-Apr;12(2):232-3.
% PMID: 4000083
%
% Constantinides CD, Atalar E, McVeigh ER.
% Signal-to-noise measurements in magnitude images from NMR phased arrays.
% Magn Reson Med. 1997 Nov;38(5):852-7.
% PMID: 9358462
%
% Kellman P, McVeigh ER.
% Image reconstruction in SNR units: a general method for SNR measurement.
% Magn Reson Med. 2005 Dec;54(6):1439-47.
% PMID: 16261576

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2009/nov/30
% Copyright © 2006-2012 Jonathan R. Polimeni and
%   The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
% $Id: mrir_array_combine_noise_correction_lookup.m,v 1.4 2009/11/30 23:19:40 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.4 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  load mrir_array_combine_noise_correction_lookup__table

  if ( isempty(mrir_array_combine_noise_correction_lookup__table{Ncha}) ),
    warning(sprintf('magnitude bias correction not yet calculated for: (( %3d channels )) -- update lookup tables', Ncha));

    snr_meas  = [];
    snr_true  = snr_true;
    std_ratio = [];
    avg_ratio = [];
    snr_bias  = [];

    return;

  end;

  snr_meas  = mrir_array_combine_noise_correction_lookup__table{Ncha}.snr_meas  ;
  snr_true  = mrir_array_combine_noise_correction_lookup__table{Ncha}.snr_true  ;
  std_ratio = mrir_array_combine_noise_correction_lookup__table{Ncha}.std_ratio ;
  avg_ratio = mrir_array_combine_noise_correction_lookup__table{Ncha}.avg_ratio ;
  snr_bias  = mrir_array_combine_noise_correction_lookup__table{Ncha}.snr_bias  ;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_combine_noise_correction_lookup.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
