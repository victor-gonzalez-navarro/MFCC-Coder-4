function [q_mfcc, Em, Xm] = mfcc_coder(Phi, win, shift_factor, x, gamma)
  % [q_mfcc, Em, Xm] = mfcc_coder(Phi, win, shift_factor, x, gamma)
  % [Computes MFCC and return quantized values
  % Input parameters:
  % - Phi: Mel filters
  % - win: window analysis
  % - shift_factor: signal is divided in frases of the window analysis length, shifted this factor of the length
  % - x: input signal (column vector)
  % - gamma: if gamma is 2 (def.), the energy is computed. But other values can be tried: sum |X(f)|^gamma
  % Return values:
  %  - q_mfcc(k,n): quantized mffc values; k coefficient of frame n 
  %  - Xm: analysis |X(f)|; Xm(k,n), |X(f)| coefficent k, frame n
  %  - Em: analysis energy band; Em(b,n), energy in band b, frame n


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

  % Split into frames and compute STFT
  xm = signal2frames(x, shift_factor, win);

%  TODO:
%  Normalize each frame with its energy
%  p = std(xm);
%  xm = xm ./ p;
%  And quantize and send p, energy of each frame

  % normalize for energy

  lenfft = size(Phi,2);
  Xm = abs(fft(xm,lenfft));
  Xm2= Xm.^gamma;
  
  % Energies in the filter bank
  Em = Phi * Xm2;

  % Missing: Transform Em into mfcc and quantize
  mfcc = dct(log(Em));

  q_mfcc = mfcc;

end
