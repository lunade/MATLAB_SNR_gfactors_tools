function [img_peft, prot] = mrir_conventional_2d(raw, varargin);
%MRIR_CONVENTIONAL_2D  reconstructs conventional (cartesian) acquisitions
%
% img_uncomb = mrir_conventional_2d(raw);
%
% img_uncomb = mrir_conventional_2d(raw, prot);
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
% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jun/01
% $Id: mrir_conventional_2d.m,v 1.2 2011/03/28 04:14:46 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  prot = [];
  DO__IMAGE_CROP = 1;
  
  if ( nargin >= 2 ), 
    prot = varargin{1}; 
    if ( ~isstruct(prot) & prot == 0 ),
      prot = [];
      DO__IMAGE_CROP = 0;
    end;
  end;
  if ( isempty(prot) ), prot = read_meas_prot__struct; end;

%  if ( mrir_ice_dimensions(raw, 'par') > 1 ),
%    warning(sprintf('input "%s" contains data along dimension "par"', inputname(1)));
%  end;


  %==--------------------------------------------------------------------==%

  hyb_roft = mrir_iDFT_freqencode(raw);
  img_peft = mrir_iDFT_phasencode(hyb_roft, 'lin', prot.lPhaseEncodingLines);

  if ( DO__IMAGE_CROP )
    img_peft = mrir_image_crop(img_peft, prot.flReadoutOSFactor);
  end


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_conventional_2d.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
