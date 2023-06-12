function snr = mrir_array_SNR_rss(sensitivity, covmtx)
%MRIR_ARRAY_SNR_RSS
%
% snr = mrir_array_SNR_rss(sensitivity, covmtx);
% snr = mrir_array_SNR_rss(img_uncombined, covmtx);

% note for "coil SNR" use sensitivities but for "image SNR" use uncombined
% images.
    
% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2007/mar/08
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
% $Id: mrir_array_SNR_rss.m,v 1.1 2007/03/23 23:41:01 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.1 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  FLAG__PREWHITEN = 0;


  %==--------------------------------------------------------------------==%

  Ncol = size(sensitivity, mrir_DIM_COL);
  Nlin = size(sensitivity, mrir_DIM_LIN);
  Ncha = size(sensitivity, mrir_DIM_CHA);
  Npar = size(sensitivity, mrir_DIM_PAR);
  Nslc = size(sensitivity, mrir_DIM_SLC);

  if ( FLAG__PREWHITEN ),
    % whitening operator for pre-whitening data
    W = mrir_array_whitening_operator(covmtx, 'svd');

    sensitivity_whitened = reshape(sensitivity, [], size(W,1)) * W.';
    sensitivity = reshape(sensitivity_whitened, size(sensitivity));

  end;

  % preallocate  1     2  3  4  5  6  7  8     9     0
  snr = zeros(Ncol, Nlin, 1, 1, 1, 1, 1, 1, Npar, Nslc);

  for ll = 1:Nslc,
    for kk = 1:Npar,
      for ii = 1:Ncol,
        for jj = 1:Nlin,
          S = squeeze(sensitivity(ii, jj, :,1,1,1,1,1,kk,ll));

          signalmag = abs( S' * S );

          if ( FLAG__PREWHITEN ),
            noisepower = abs(S' * S);
          else,
            noisepower = abs(S' * covmtx * S);
          end;

          snr(ii, jj, 1,1,1,1,1,1,kk,ll) = signalmag / sqrt(noisepower);
        end;
      end;
    end;
  end;

  % if noise covariance is calculated in k-space a scale factor may need to
  % be applied to account for averaging in FFT. depends on convention used.
  % typically it is best to have a "image-domain" noise covariance that has
  % already been scaled to account for the FFT.
  scale_factor = 1 / sqrt(2 * Ncol * Nlin * Npar);

  snr = scale_factor * snr;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/COIL_ARRAYS/MATLAB/mrir_array_SNR_rss.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
