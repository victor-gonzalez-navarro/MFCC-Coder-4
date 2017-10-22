function [xr D Emq Xmr] = mfcc_decoder(PhiI, win, shift_factor, mfcc, max_ite, D_threshold, intermediate_frames, gamma)
  % [xr D Emq Xmr] = mfcc_decoder(PhiI, win, shift_factor, mfcc, max_ite, D_threshold, intermediate_frames, gamma)
  % Decoded the signal xr from the quantized mffc values
  % Input parameters:
  % - PhiI: Pseudo inverse of mel filters
  % - win: window analysis
  % - shift_factor: signal is divided in frases of the window analysis length, shifted this factor of the length
  % - mfcc: quantized mffc values; k coefficient of frame n 
  % - max_ite (def. 100), D_threshold (def. 1e-6), control the convergence of iterative algorithm lsee_mstft
  % - intermediate_frames: creates interpolated frames, in energy domain, to help lsee_mstft  
  % - gamma: if gamma is 2 (def.), the energy is computed. But other values can be tried: sum |X(f)|^gamma
  % Return values:
  %  - xr: decoded signal 
  %  - D: vector with the distorsion at each iteration (to analyze convergence)
  %  - Emq: Dequantized energy band; Emq(b,n), energy in band b, frame n
  %  - Xmr: Reconstructed |X(f)|; Xm(k,n), |X(f)| coefficent k, frame n


  % ---------------------------------------------------------------------
  % Copyright (C) Antonio Bonafonte, 2016
  % Universitat Politecnica de Catalunya, Barcelona, Spain.
  % 
  % Permission to copy, use, modify, sell and distribute this software
  % is granted provided this copyright notice appears in all copies.
  % This software is provided "as is" without express or implied
  % warranty, and with no claim as to its suitability for any purpose.
  % 
  % ---------------------------------------------------------------------



  % TO DO: q^-1 
  % For the moment, it arrives q(mfcc)

  Emq = exp(idct(mfcc));

  if nargin < 5
     max_ite = 100;
  end
  if nargin < 6
     D_threshold = 1e-4;
  end
  if nargin < 7
     intermediate_frames = 0;
  end
  if nargin < 8
     gamma = 2;
  end

  % Interpolate energies in intermediate frames
  if intermediate_frames > 0
    nframes = size(Emq,2);
    idx=1:nframes;
    step = 1/(intermediate_frames+1);
    idxn = 1:step:nframes;

    Etmp = interp1(idx, Emq', idxn);
    Emqi = Etmp';
  else
    Emqi = Emq;
  end

  % Envelope approximation
  Xmr2 = PhiI * Emqi;
  Xmr  = Xmr2 .^ (1/gamma);

  % Estimate signal x, that match Xmr
  lenfft = size(PhiI, 1);
  [xr D] = lsee_mstft(Xmr, win, win, shift_factor/(intermediate_frames+1), lenfft, max_ite, D_threshold);
end
