function ft = mrir_iDFT(raw, dim, varargin)
%MRIR_IDFT  inverse Discrete Fourier Transform
%
% ft = mrir_iDFT(raw, dim)
% ft = mrir_iDFT(raw, dim, N)
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
% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jan/24
% $Id: mrir_iDFT.m,v 1.4 2011/03/28 04:14:46 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.4 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  if ( nargin >= 3 ),
    Npoint = varargin{1};
  else,
    Npoint = size(raw, dim);
  end;


  %==--------------------------------------------------------------------==%

  if ( Npoint > size(raw, dim) ),

    % because k-space ordering requires FFT shifting, we must zero pad
    % before the shift, thus simply computing an N-point FFT will not
    % suffice.
    pad_dims = zeros(1, 16);

    pad_dims(dim) = Npoint - size(raw, dim);

    % by convention, always pad AFTER the data
    raw = mrir_zeropad(raw, pad_dims, 'post');

  end;

  ft = fftshift(ifft(ifftshift(raw, dim), Npoint, dim), dim) * Npoint;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_iDFT.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:




  % note on scaling:
  %
  %  in here and in "mrir_fDFT.m", i have reversed the MATLAB convention for
  %  scaling the fft. MATLAB does not like the 'unitary' normalization that
  %  places the \sqrt{N} before both the fDFT and iDFT. see the help for
  %  "fft.m":
  %
  %%       For length N input vector x, the DFT is a length N vector X,
  %%       with elements
  %%                        N
  %%          X(k) =       sum  x(n)*exp(-j*2*pi*(k-1)*(n-1)/N), 1 <= k <= N.
  %%                       n=1
  %%       The inverse DFT (computed by IFFT) is given by
  %%                        N
  %%          x(n) = (1/N) sum  X(k)*exp( j*2*pi*(k-1)*(n-1)/N), 1 <= n <= N
  %%                       k=1
  %
  %  it would have been better if i originally stuck with the unitary
  %  version, but other code (e.g., the SNR calculation code) assumes the
  %  scaling convention established in this function, so i am afraid to
  %  change it now. i have given the 1/N to the forward, and 1 to the
  %  inverse, so that still the composition of the two yields the identity,
  %  but this makes comparison with operators using the unitary
  %  scaling---like the DFT matrix.
