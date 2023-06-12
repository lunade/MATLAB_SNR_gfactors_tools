function img_combine_rss = mrir_array_combine_rss(img_multichan, varargin)
%MRIR_ARRAY_COMBINE_RSS
%
% img_combine_rss = mrir_array_combine_rss(img_multichan);

% TODO: return measure of phase compatible with root-sum-of-squares
% combination, e.g., homodyne detection (Noll et al., 1991)
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
% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jan/27
% $Id: mrir_array_combine_rss.m,v 1.2 2011/03/28 04:14:45 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.2 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

%   alpha = 1;
%
%   if ( nargin >= 2 ),
%     covmtx = varargin{1};
%
%     % arbitrary constant (normalizes image to approx same level as RSS combo)
%     alpha = sqrt(mean(diag(covmtx)));
%
%     W = mrir_array_whitening_operator(covmtx, 'svd');
%     img = mrir_array_whitening_apply(img_multichan, W);
%
%   else,
%     img = img_multichan;
%   end;

%  img_combine_rss = alpha * squeeze(sqrt(sum(abs(img).^2,3)));
%  img_combine_rss = alpha * (sqrt(sum(abs(img).^2,3)));
  img_combine_rss = (sqrt(sum(abs(img_multichan).^2,3)));


  %  img_combine_phz = angle(sum(img_multichan, 3));

%  img_combine_cplx = img_combine_rss .* exp(i * img_combine_phz);

%  if ( nargout > 1 ),

%    varargout{1} = img_combine_phz;
%    varargout{2} = img_combine_cplx;

%  end;



  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_combine_rss.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End: