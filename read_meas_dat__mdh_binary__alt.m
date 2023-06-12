function [mdh, scan_num] = read_meas_dat__mdh_binary(binary, varargin)
%READ_MEAS_DAT__MDH_BINARY  quickly parse MDH extracted as one binary vector
%
% mdh_struct = read_meas_dat__mdh_binary(binary)

% Copyright © 2006-2012 Jonathan R. Polimeni and
%   The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense


% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2008/mar/24
% $Id: read_meas_dat__mdh_binary__alt.m,v 1.7 2012/05/26 09:58:27 jonp Exp $
%
% VD-platform support:
% Anastasia Yendiki <ayendiki@nmr.mgh.harvard.edu>, 2011/oct/24
%**************************************************************************%

  VERSION = '$Revision: 1.7 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  % allow user to request full MDH header or an abbreviated version
  request_full = 0;
  if ( nargin >= 2 ),
    request_full = varargin{1};
  end;

  IS__VD_PLATFORM = 0;
  if ( length(binary) == 192 ),  % Is this a VD11 data file?
    IS__VD_PLATFORM = 1;
  end;


  if ( IS__VD_PLATFORM ),                      % VD11 or newer MDH structure
    %==------------------------------------------------------------------==%

    % constants defined in <n4/pkg/MrServers/MrMeasSrv/SeqIF/MDH/mdh.h>
    MDH_NUMBEROFEVALINFOMASK   = 2;
    MDH_NUMBEROFICEPROGRAMPARA = 24;  % value differs from VB-platform

    MDH_FREEHDRPARA = 4;


    %==------------------------------------------------------------------==%

    ind = 1;

    mdh.ulDMALength = typecast(binary(ind:ind+1), 'uint16');
    ind = ind + 2;

    ulFlags1    = typecast(binary(ind:ind+0), 'uint8');
    ind = ind + 1;
    ulFlags2    = typecast(binary(ind:ind+0), 'uint8');
    ind = ind + 1;

    mdh.ulFlagsAndDMALength        = double(ulFlags2) * 2^24 + ...
                                     uint32(double(ulFlags1) * 2^16) + ...
                                     uint32(mdh.ulDMALength);

    mdh.lMeasUID                   = typecast(binary(ind:ind+3), 'int32');
    ind = ind + 4;
    mdh.ulScanCounter              = typecast(binary(ind:ind+3), 'uint32');
    ind = ind + 4;
    mdh.ulTimeStamp                = typecast(binary(ind:ind+3), 'uint32'); % 2.5 ms ticks since 00:00
    mdh.flTimeStamp_ms             = double(mdh.ulTimeStamp) * 2.5;
    ind = ind + 4; % milliseconds
    mdh.ulPMUTimeStamp             = typecast(binary(ind:ind+3), 'int32');  % 2.5 ms ticks since last trigger
    mdh.flPMUTimeStamp_ms          = double(mdh.ulPMUTimeStamp) * 2.5;
    ind = ind + 4; % milliseconds
    mdh.ushSystemType              = typecast(binary(ind:ind+1), 'uint16');
    ind = ind + 2;

    if ( request_full ),
      mdh.ushPTABPosDelay          = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.lPTABPosX                = typecast(binary(ind:ind+3), 'int32');
      ind = ind + 4;
      mdh.lPTABPosY                = typecast(binary(ind:ind+3), 'int32');
      ind = ind + 4;
      mdh.lPTABPosZ                = typecast(binary(ind:ind+3), 'int32');
      ind = ind + 4;
      mdh.ulReserved1              = typecast(binary(ind:ind+3), 'uint32');
      ind = ind + 4;
    else
      ind = ind + 18;
    end;

    mdh.aulEvalInfoMask(1:MDH_NUMBEROFEVALINFOMASK) = typecast(binary(ind:ind+4*MDH_NUMBEROFEVALINFOMASK-1), 'uint32');
    ind = ind + 4*MDH_NUMBEROFEVALINFOMASK;

    % build 64-bit mask from two 32-bit integers
    mask = double(...
        bitshift(uint64(mdh.aulEvalInfoMask(2)), 32)) ...
                  + double(mdh.aulEvalInfoMask(1));

    mdh.ushSamplesInScan           = typecast(binary(ind:ind+1), 'uint16');
    ind = ind + 2;
    mdh.ushUsedChannels            = typecast(binary(ind:ind+1), 'uint16');
    ind = ind + 2;

    if (1),
      %  sLoopCounter
      mdh.ushLine                  = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushAcquisition           = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushSlice                 = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushPartition             = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushEcho                  = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushPhase                 = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushRepetition            = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushSet                   = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushSeg                   = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushIda                   = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushIdb                   = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushIdc                   = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushIdd                   = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushIde                   = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
    end;

    if (1),
      % sCutOffData
      mdh.ushPre                   = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushPost                  = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
    end;

    mdh.ushKSpaceCentreColumn      = typecast(binary(ind:ind+1), 'uint16');
    ind = ind + 2;
    mdh.ushCoilSelect              = typecast(binary(ind:ind+1), 'uint16');
    ind = ind + 2;
    mdh.fReadOutOffcentre          = typecast(binary(ind:ind+3), 'single');
    ind = ind + 4;
    mdh.ulTimeSinceLastRF          = typecast(binary(ind:ind+3), 'uint32');
    ind = ind + 4;
    mdh.ushKSpaceCentreLineNo      = typecast(binary(ind:ind+1), 'uint16');
    ind = ind + 2;
    mdh.ushKSpaceCentrePartitionNo = typecast(binary(ind:ind+1), 'uint16');
    ind = ind + 2;

    if (1),
      % sSliceData
      if (1),
        % sVector
        mdh.flSag                  = typecast(binary(ind:ind+3), 'single');
        ind = ind + 4;
        mdh.flCor                  = typecast(binary(ind:ind+3), 'single');
        ind = ind + 4;
        mdh.flTra                  = typecast(binary(ind:ind+3), 'single');
        ind = ind + 4;
      end;
      if ( request_full ),
        mdh.aflQuaternion(1:4)     = typecast(binary(ind:ind+15), 'single');
      end;
      ind = ind + 16;
    end;

    if ( request_full ),
      mdh.aushIceProgramPara(1:MDH_NUMBEROFICEPROGRAMPARA) = ...
        typecast(binary(ind:ind+2*MDH_NUMBEROFICEPROGRAMPARA-1), 'uint16');
      ind = ind + 2*MDH_NUMBEROFICEPROGRAMPARA;

      mdh.aushFreePara(1:MDH_FREEHDRPARA) = ...
        typecast(binary(ind:ind+2*MDH_FREEHDRPARA-1), 'uint16');
      ind = ind + 2*MDH_FREEHDRPARA;
    else,
      ind = ind + 2*(MDH_NUMBEROFICEPROGRAMPARA + MDH_FREEHDRPARA);
    end;

    mdh.ushApplicationCounter      = typecast(binary(ind:ind+1), 'uint16');
    ind = ind + 2;
    mdh.ushApplicationMask         = typecast(binary(ind:ind+1), 'uint16');
    ind = ind + 2;
    mdh.ulCRC                      = typecast(binary(ind:ind+3), 'uint32');
    ind = ind + 4;

    mdh.ushChannelId               = [];       % for backwards compatibility

    %======================================================================%
    %======================================================================%
  else,                                          % VB-platform MDH structure
    %======================================================================%
    %======================================================================%

    % constants defined in <n4/pkg/MrServers/MrMeasSrv/SeqIF/MDH/mdh.h>
    MDH_NUMBEROFEVALINFOMASK   = 2;
    MDH_NUMBEROFICEPROGRAMPARA = 4;  % value differs from VD-platform

    MDH_FREEHDRPARA = 4;


    %==------------------------------------------------------------------==%

    ind = 1;

    mdh.ulDMALength = typecast(binary(ind:ind+1), 'uint16');
    ind = ind + 2;

    ulFlags1    = typecast(binary(ind:ind+0), 'uint8');
    ind = ind + 1;
    ulFlags2    = typecast(binary(ind:ind+0), 'uint8');
    ind = ind + 1;

    mdh.ulFlagsAndDMALength        = double(ulFlags2) * 2^24 + ...
                                     uint32(double(ulFlags1) * 2^16) + ...
                                     uint32(mdh.ulDMALength);

    mdh.lMeasUID                   = typecast(binary(ind:ind+3), 'int32');
    ind = ind + 4;
    mdh.ulScanCounter              = typecast(binary(ind:ind+3), 'uint32');
    ind = ind + 4;
    mdh.ulTimeStamp                = typecast(binary(ind:ind+3), 'uint32'); % 2.5 ms ticks since 00:00
    mdh.flTimeStamp_ms             = double(mdh.ulTimeStamp) * 2.5;
    ind = ind + 4; % milliseconds
    mdh.ulPMUTimeStamp             = typecast(binary(ind:ind+3), 'int32');  % 2.5 ms ticks since last trigger
    mdh.flPMUTimeStamp_ms          = double(mdh.ulPMUTimeStamp) * 2.5;
    ind = ind + 4; % milliseconds

    mdh.aulEvalInfoMask(1:MDH_NUMBEROFEVALINFOMASK) = typecast(binary(ind:ind+4*MDH_NUMBEROFEVALINFOMASK-1), 'uint32');
    ind = ind + 4*MDH_NUMBEROFEVALINFOMASK;

    % build 64-bit mask from two 32-bit integers
    mask = double(...
        bitshift(uint64(mdh.aulEvalInfoMask(2)), 32)) ...
                  + double(mdh.aulEvalInfoMask(1));

    mdh.ushSamplesInScan           = typecast(binary(ind:ind+1), 'uint16');
    ind = ind + 2;
    mdh.ushUsedChannels            = typecast(binary(ind:ind+1), 'uint16');
    ind = ind + 2;

    if (1),
      %  sLoopCounter
      mdh.ushLine                  = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushAcquisition           = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushSlice                 = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushPartition             = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushEcho                  = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushPhase                 = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushRepetition            = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushSet                   = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushSeg                   = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushIda                   = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushIdb                   = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushIdc                   = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushIdd                   = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushIde                   = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
    end;

    if (1),
      % sCutOffData
      mdh.ushPre                   = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
      mdh.ushPost                  = typecast(binary(ind:ind+1), 'uint16');
      ind = ind + 2;
    end;

    mdh.ushKSpaceCentreColumn      = typecast(binary(ind:ind+1), 'uint16');
    ind = ind + 2;
    mdh.ushDummy                   = typecast(binary(ind:ind+1), 'uint16');
    ind = ind + 2;
    mdh.fReadOutOffcentre          = typecast(binary(ind:ind+3), 'single');
    ind = ind + 4;
    mdh.ulTimeSinceLastRF          = typecast(binary(ind:ind+3), 'uint32');
    ind = ind + 4;
    mdh.ushKSpaceCentreLineNo      = typecast(binary(ind:ind+1), 'uint16');
    ind = ind + 2;
    mdh.ushKSpaceCentrePartitionNo = typecast(binary(ind:ind+1), 'uint16');
    ind = ind + 2;

    if ( request_full ),
      mdh.aushIceProgramPara(1:MDH_NUMBEROFICEPROGRAMPARA) = ...
        typecast(binary(ind:ind+2*MDH_NUMBEROFICEPROGRAMPARA-1), 'uint16');
      ind = ind + 2*MDH_NUMBEROFICEPROGRAMPARA;

      mdh.aushFreePara(1:MDH_FREEHDRPARA) = ...
        typecast(binary(ind:ind+2*MDH_FREEHDRPARA-1), 'uint16');
      ind = ind + 2*MDH_FREEHDRPARA;
    else
      ind = ind + 2*(MDH_NUMBEROFICEPROGRAMPARA + MDH_FREEHDRPARA);
    end;

    if (1),
      % sSliceData
      if (1),
        % sVector
        mdh.flSag                  = typecast(binary(ind:ind+3), 'single');
        ind = ind + 4;
        mdh.flCor                  = typecast(binary(ind:ind+3), 'single');
        ind = ind + 4;
        mdh.flTra                  = typecast(binary(ind:ind+3), 'single');
        ind = ind + 4;
      end;
      if ( request_full ),
        mdh.aflQuaternion(1:4)     = typecast(binary(ind:ind+15), 'single');
      end;
      ind = ind + 16;
    end;

    mdh.ushChannelId               = typecast(binary(ind:ind+1), 'uint16');
    ind = ind+2;

    if ( request_full ),
      mdh.ushPTABPosNeg            = typecast(binary(ind:ind+1), 'uint16');
    end;
    ind = ind+2;

    %==------------------------------------------------------------------==%
  end;

  scan_num = double(mdh.ulScanCounter);

  return;



  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/read_meas_dat__mdh_binary__alt.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:

