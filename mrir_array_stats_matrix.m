function varargout = mrir_array_stats_matrix(rawdata, varargin)
%MRIR_ARRAY_STATS_MATRIX
%
% statmtx = mrir_array_stats_matrix(rawdata, stat_type)
%
%  'cor':  correlation matrix
%  'cof':  correlation coefficient matrix
%  'cov':  covariance matrix
%  'std':  standard deviation matrix
%  'avg':  average vector
%
%
% rawdata must be stored such that the THIRD dimension corresponds to the
% coil channels, but the other dimensions can be ordered arbitrarily.
%
% see also MRIR_ARRAY_STATS_MATRIX__PERVOXEL.


% statmtx = mrir_array_stats_matrix(rawdata, stat_type, do_bandwidth_scale)

% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/mar/08
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
% $Id: mrir_array_stats_matrix.m,v 1.7 2011/08/20 22:16:06 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.7 $';
  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  if ( nargin >= 2 ),
    stat_type = varargin{1};
  else,
    stat_type = 'cov';  % covariance calculation by default
  end;


  if ( nargin >= 3 ),
    DO__NOISE_BW_SCALING = varargin{2};
  else,
    DO__NOISE_BW_SCALING = 0;
  end;


  % assume here that raw data contains noise, which should be arranged as
  % a collection of lines across all channels

  % TODO: generalize this to work with signal stats matrices, where there
  % may be other dimensions (e.g., slices) in the data

  % raw data from scanner is typically stored as single-precision floats,
  % but double-precision is better for calculations.

  if ( isa(rawdata, 'single') ),
%    disp('converting single-precision input data to double-precision floats for processing');
  end;

  noise = double(rawdata);

  if ( DO__NOISE_BW_SCALING ),
    noise_bandwidth = mrir_noise_bandwidth(noise);

    if ( noise_bandwidth < 0.6 ),
      warning('noise bandwidth is too low; this data is unlikely to be pure noise');
    end;

  end;

  dims = size(noise);
  if ( length(dims) == 2 ),
    dims = [dims(1) dims(2), 1];
  end;


  % by convention, ICE data is stored as cols x rows x chan x ..., so reshape to
  % samples X chan
  order = 1:length(dims);
  order(3) = length(dims);
  order(end) = 3;

  % move channel index to end
  noise = permute(noise, order);

  noise = reshape( noise, [], dims(3) );

  % transpose so each observation is a column vector
  noise = permute(noise, [2, 1]);

  if ( nnz(noise) ~= numel(noise) ),
    warning(sprintf(' %2.1f%% of noise data elements are zero-valued; statistics may be inaccurate', [1 - (nnz(noise)/numel(noise))]*100));
  end;

  Nchan = size(noise, 1);
  Nsamp = size(noise, 2);

  switch lower(stat_type),
   case 'mtx',
    varargout{1} = noise;
   case 'cor',
    % normalize the correlation matrix such that it is equivalent to
    % covariance matrix in the case of truly zero-mean noise
    cormtx = ( noise*noise' )                          / (Nsamp-1);

    if ( DO__NOISE_BW_SCALING ),
      disp(sprintf('\nscaling covariance matrix by effective noise bandwidth:  %.4f\n', noise_bandwidth));

      % scale by noise bandwidth (Kellman & McVeigh, 2005)
      cormtx = cormtx / noise_bandwidth;
    end;

    varargout{1} = cormtx;
   case 'cov',
    avgvec = mean(noise, 2);
    covmtx = ( noise*noise' - (avgvec*avgvec')*Nsamp ) / (Nsamp-1);
    %covmtx = cov(noise.').';

    if ( DO__NOISE_BW_SCALING ),
      disp(sprintf('\nscaling covariance matrix by effective noise bandwidth:  %.4f\n', noise_bandwidth));

      % scale by noise bandwidth (Kellman & McVeigh, 2005)
      covmtx = covmtx / noise_bandwidth;
    end;

    varargout{1} = covmtx;
   case 'cof',
    avgvec = mean(noise, 2);
    covmtx = ( noise*noise' - (avgvec*avgvec')*Nsamp ) / (Nsamp-1);
    varvec = diag(covmtx);
    varmtx = diag(varvec);
    stdmtx = sqrt(varmtx);
    rotmtx = inv(stdmtx);
    cofmtx = rotmtx * covmtx * rotmtx;
    %cofmtx = corrcoef(noise.').';
    varargout{1} = cofmtx;
   case 'std',
    avgvec = mean(noise, 2);
    covmtx = ( noise*noise' - (avgvec*avgvec')*Nsamp ) / (Nsamp-1);
    varvec = diag(covmtx);
    stdvec = sqrt(varvec);
    stdsum = sum(stdvec);
    stdvec = stdvec / stdsum * Nchan;
    varargout{1} = stdvec;
   case {'mean', 'avg'},
    avgvec = mean(noise, 2);
    varargout{1} = avgvec;
   case 'avgmtx',
    avgvec = mean(noise, 2);
    avgmtx = avgvec*avgvec';
    varargout{1} = avgmtx;
   case 'res',
    % compute resolution matrix (like point-spread, should be identity if pinv is true inv)
    resmtx = ( noise*pinv(noise) )                     / (Nsamp-1);

    varargout{1} = resmtx;
   case 'onc',

    % ortho-normalized correlation (Bruehrer et al. 2007)
    oncmtx = 0;
    for ind = 1:size(noise,2),
      oncmtx = oncmtx + ( noise(:,ind) * pinv(noise(:,ind)) );
    end;

    varargout{1} = oncmtx;

   case 'snr',

    % recursive!
%    noisecov = mrir_array_stats_matrix(rawdata, 'cov', 1);
%
%    W = mrir_array_whitening_operator(noisecov, 'svd');


   otherwise,
    error('unknown statistics matrix type: %s', stat_type);
  end;



  if ( nargout >= 2 ),
    varargout{2} = noise;
  end;



  % mean(abs(covmtx(find(tril(ones(size(covmtx)), -1)))))


%%   % alternative equivalent formulae:
%%
%%   %% Pruessmann's recipe for COVARIANCE, corrected version:
%%
%%   for jj = 1:Nchan,
%%     eta_j = noise(:,jj);
%%     for kk = 1:Nchan,
%%       eta_k = noise(:,kk);
%%       prumtx_new(jj,kk) = 0.5 * [var(eta_j + eta_k) + i*var(eta_j + i*eta_k) ...
%%                    - (1+i)*( var(eta_j) + var(eta_k) )];
%%     end;
%%   end;
%%
%%
%%   %% matrix operations for COVARIANCE:
%%
%%   m = mean(Z);
%%   M = repmat(m, Nsamp, 1);
%%
%%   covmtx = ( Z'*Z + M'*M - Z'*M - M'*Z ) / (Nsamp-1);
%%   covmtx = ( Z'*Z - M'*M )               / (Nsamp-1);
%%   covmtx = ( Z'*Z - (m'*m)*Nsamp )       / (Nsamp-1);
%%
%%
%%   %% matrix operations for CORRELATION COEFFICIENT:
%%   vars = diag(covmtx);
%%   varm = diag(vars);
%%   s = sqrt(varm);
%%   w = inv(s);
%%   cofmtx = w*covmtx*w;




  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/mrir_array_stats_matrix.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
