function roft = mrir_iDFT_freqencode(raw, varargin)
%MRIR_IDFT_FREQENCODE  inverse Discrete Fourier Transform along freq encode
%
% roft = mrir_iDFT_freqencode(raw)
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
% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jan/30
% $Id: mrir_iDFT_freqencode.m,v 1.2 2011/03/28 04:14:46 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  Npoint = size(raw, 1);
  if ( nargin >= 2 && ~isempty(varargin{2-1}) ),
    Npoint = varargin{2-1};
  end;

  roft = mrir_iDFT(raw, 1, Npoint);

  
  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_iDFT_freqencode.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
