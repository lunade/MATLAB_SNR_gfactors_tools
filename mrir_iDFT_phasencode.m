function peft = mrir_iDFT_phasencode(raw, coordinate_str, varargin)
%MRIR_IDFT_PHASENCODE  inverse Discrete Fourier Transform along phase encode
%
% peft = mrir_iDFT_phasencode(raw, coordinate_str)
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
% $Id: mrir_iDFT_phasencode.m,v 1.2 2011/03/28 04:14:46 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  % assume primary phase encode direction by default
  if ( ~exist('coordinate_str', 'var') ),
    coordinate_str = 'lin';
  end;

  switch lower(coordinate_str(1:3)),
   case 'lin',
    dim = 2;
   case 'par',
    dim = 9;
   otherwise,
    error('unrecognized data dimension: "%s"', coordinate_str);
  end;

  if ( size(raw, dim) == 1 ),
    warning(sprintf('input "%s" contains no data along dimension "%s"', ...
		    inputname(1), coordinate_str));
  end;

  Npoint = size(raw, dim);
  if ( nargin >= 3 && ~isempty(varargin{3-2}) ),
    Npoint = varargin{3-2};
  end; 
  
  peft = mrir_iDFT(raw, dim, Npoint);
  

  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_iDFT_phasencode.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
  