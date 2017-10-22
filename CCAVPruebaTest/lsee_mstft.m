function [x D] = lsee_mstft(Ym, winA, winS, shift_factor, nfft, max_ite, threshold)
  % [x D] = lsee-mstft(Xnr, winA, winS, shift_factor, max_ite, D_threshold)
  % Compute signal x(n) so that the magnitude of its STFT is as similar to the magnitude given in Ym(freq, frame)
  %
  % Griffin and Lim, Signal Estimation from Modified STFT, ASSP-32, 2, 1984


  % ---------------------------------------------------------------------
  % Copyright (C) Antonio Bonafonte, 2014
  % Universitat Politecnica de Catalunya, Barcelona, Spain.
  % 
  % Permission to copy, use, modify, sell and distribute this software
  % is granted provided this copyright notice appears in all copies.
  % This software is provided "as is" without express or implied
  % warranty, and with no claim as to its suitability for any purpose.
  % 
  % ---------------------------------------------------------------------

  len = length(winA);
  D   = zeros(1,max_ite);

%  Dh  = zeros(1,max_ite);
%  R   = zeros(1,max_ite);

  winNorm = winA .* winS;
  xmr = real(ifft(Ym));   % Should be real as Ym is real (|.|)
  xmr = xmr(1:len,:);     % from nftt to len

  z = frames2signal(xmr, shift_factor, winS, 0, winNorm);

  %xr = randn(size(z));  %Init aleatoria
  xr = z;
  for i=1:max_ite
    xmr = signal2frames(xr, shift_factor, winA);
    Xmr = fft(xmr,nfft);
    D(i) = sum(sum((abs(Xmr) - abs(Ym)).^2))/numel(Ym);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Zm = Ym .* exp(j * angle(Xmr));

% Debug ....
%    zmC = ifft(Zm);
%    zmCr = real(zmC); sumR = sum(abs(zmCr(:)));
%    zmCi = imag(zmC); sumI = sum(abs(zmCi(:)));
%    R(i) = sumI/sumR;
%    Dh(i) = sum(sum((abs(Xmr(1:nfft/2,:)) - abs(Ym(1:nfft/2,:))).** 2))/(numel(Ym)/2);

    zm = real(ifft(Zm));
    zm = zm(1:len,:);
    xr = frames2signal(zm, shift_factor, winS, 0, winNorm);

    if (i > 1 && i < max_ite && (D(i-1)-D(i))/D(i-1) < threshold)
      D(i+1:end) = [];
%      R(i+1:end) = [];
      break
    end
%    fprintf('D(%d) = %f\tDh(%d) = %f\tI/R(%d)= %f\n', i, D(i), i, Dh(i), i, R(i)); fflush(stdout);
    %fprintf('D(%d) = %f\n', i, D(i)); fflush(stdout);?????????????????????????????????????????????????????????????????????????????????????????????????
  end
  x = xr;
end
