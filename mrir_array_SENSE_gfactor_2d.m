function gfactor = mrir_array_SENSE_gfactor_2d(sensitivity, covmtx, R1, R2)
%MRIR_ARRAY_SENSE_GFACTOR_2D  estimates SENSE g-factor for two acceleration directions
%
% gfactor = mrir_array_SENSE_gfactor_2d(sensitivity, covmtx, R1, R2)
%
%
% see also MRIR_ARRAY_SENSE_GFACTOR_2D__REGULARIZE.
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

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/apr/04
% $Id: mrir_array_SENSE_gfactor_2d.m,v 1.5 2012/03/03 21:30:55 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.5 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  Ncol = size(sensitivity, 1); % frequency encoded
  Nlin = size(sensitivity, 2); % phase encoded
  Ncha = size(sensitivity, 3);

  covmtxinv = inv(covmtx);

  alias_map1 = mrir_array_accelerated_aliasmap(Nlin, R1);
  alias_map2 = mrir_array_accelerated_aliasmap(Ncol, R2);

%  if ( prod(size(alias_map1)) ~= Nlin ),
%    warning('R1 not compatible with Nlin');
%  end;

%  if ( prod(size(alias_map2)) ~= Ncol ),
%    warning('R2 not compatible with Ncol');
%  end;

  % preallocate
  gfactor  = zeros(floor(Ncol/R2)*R2, floor(Nlin/R1)*R1);

  for jj = 1:floor(Nlin/R1),

    % indices of the aliased pixels for this position from lookup table
    aliased_ind1 = alias_map1(jj, 1:R1);

    for ii = 1:floor(Ncol/R2),

      % indices of the aliased pixels for this position from lookup table
      aliased_ind2 = alias_map2(ii, 1:R2);

      % R2 x R1 x Nchan
      b = sensitivity(aliased_ind2, aliased_ind1, :);

      % encoding matrix: Nchan x R1*R2, or Nchan x Nalias
      E = permute(reshape(b, [R1*R2, Ncha]), [2, 1]);

      nu = E' * covmtxinv * E;
      eta = inv(nu);

      gfactor(aliased_ind2,aliased_ind1) = reshape(sqrt(abs(diag(nu).*diag(eta))), [R2,R1]);

    end;
  end;


  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_SENSE_gfactor_2d.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
