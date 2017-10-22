
% Implementation based on the Paper:
%
% Low Bit-Rate Speech Coding Through Quantization of Mel-Frequency Cepstral Coefficients
% Laura E. Boucheron, Member, IEEE, Phillip L. De Leon, Senior Member, IEEE, and Steven Sandoval

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

% ALGORITHM:

%  CODER:
%  Foreach (overlapped) frames
%  1. Get frames (xm)
%  2. Compute fft: X
%  3. Compute |X| X and |X|^2  (Xm2)
%  4. Compute Energy in each mel-frequency band (Em)
%  5. Compute MFCC  ~~~~> to be done
%  6. Quantize MFCC ~~~~> to be done

% DECODER
%  1. Dequantize MFCC                      ~~~~~~> to be done
%  2. Transform MFCC to Energy bands (Emq) ~~~~~~> to be done
%  3. Interpolate Energy bands (to experiment)  (Emq)
%  4. Transform energy bands (Emq) to estimate envelope (Xmr2 and Xmr)
%  5. Tranform envelope Xmr into signal (xr[n], and Xr(f))

% ---------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration and definition of mel filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename
%basefile='ona8cs';
basefile='SA000S17';


show_frame = 100; % show this frame (<= 0: do not show)
show_mel   = 40;  % show this filter <= 0: do not show
save_plots = 1;

% Analysis options
fs       = 8000;                 % Sampling frequency
frameDur = 30e-3;                % Frame duration
len      = round(frameDur * fs); % Frame len, in samples
lenfft   = 2^ceil(log2(len));    % FFT size, ej, len=240, lenfft=256
win      = hamming(len);         % Hamming Window
gamma    = 2;                    % gamma=2 => energy band; but in ASR, gamma=1 or gamma=1.3 gives better results.

% Shift of frames.
% 1        => frames do not overlap, few parameters to quantize ...
%             but horrible quality (phase cannot be estimated)
% aprox. 0 => frames overlap a lot, 
%             better quality but many parameters to be sent
shift_factor = 0.5;

% interpolate_nframes: in the decoder, interpolate this number between each Enery bands, to have more frames.
% This may help to the LSE-ISTFTM algorithm.
% Maybe it is possible higher shift_factor and interpolation
interpolate_nframes=0; 

% Filter bank options
nfilters    = 60;%60
nfilters1k  = 30;%30

% lsee_mstftm options
max_ite     = 100;
D_threshold = 1e-6;


% Computer mel filters: each row (1..nfilters) is one filter (defined in freq.)
% Length: nftt
% Given x => |X| => Phi*|X| is a vector with the energy in each 'nfilter' bands.
[Phi, fc] = melfilters(lenfft, fs, nfilters, nfilters1k, 1);

% Compute pseudo inverse: 
PhiI  = pinv(Phi);

% Short review of pseudo inverse:
% Dimension of Phi is: nfilters x nfft. It does not have inverse.
% There are many |X| that produce the same energy bands.
% The pseudo-inverse is the best we can do: we get |X| with same
% energy bands: from all the options, the one with min. norm is selected.
% E.g.
% Given a vector x, and the matrix A, test:
% x=[1:3]', A=rand(2,3); y=A*x; xr = pinv(A)*y,  y,  yr=A*xr, norm_x=norm(x), norm_xr=norm(xr)
% xr is not x, but A*x is the same than A*xr. And the norm of xr is the smallest from all vectors with same A*x



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S I G N A L: (this can be loop for any number of files)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileIn= [basefile, '.wav']
fileOut=[basefile, 'out.wav'];

% Read signal
[x, fm] = audioread(fileIn);
if fm ~= fs
  disp 'Error: sampling freq. does not match configuration; either change configuration parameter or change rate of file'
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C O D E R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[q_mfcc Em Xm] = mfcc_coder(Phi, win, shift_factor, x, gamma);
% The coded signal is just q_mfcc. Em and Xm are just for debug and print results. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q U A N T I Z E R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[q_mfcc] = mfcc_quantizer(q_mfcc, fmatrix);
%disp('He sacado el quantizer!');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D E C O D E R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xr D Emq Xmr] = mfcc_decoder(PhiI, win, shift_factor, q_mfcc, max_ite, D_threshold, interpolate_nframes, gamma);
% The decoded signal is xr. D, Xmr and Emq are just for monitorization of the algorithms



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%audioread(xr, fs, fileOut);%??????????????????????????????????????????????????????????????????????????????????????????????
%myplay(x, fm)
sound(x,8000);
pause(3);
sound(xr,8000);

close all
set (0, 'defaultlinelinewidth', 1) % line width in plots

%error = x-xr(1:(size(x,1)),1);

t = (0:length(x)-1)/fs;
subplot(2,1,1), plot(t,x); title 'Original x(t)'
t = (0:length(xr)-1)/fs;
subplot(2,1,2), plot(t,xr), title 'Reconstructed xr(t)', 
%t = (0:length(error)-1)/fs;
%subplot(3,1,3), plot(t,error); title 'Error er(t)', 
pause(100);

%if save_plots, print  1_x_and_xr.pdf, end
subplot(1,1,1)

f       =(0:lenfft-1) * fs/lenfft;
nfft2  = lenfft/2; 
f2     = f(1:nfft2);

if show_mel > 0
  % Plot Phi
  subplot(2,1,1)
  plot(f2, Phi(:,1:nfft2));
  title(sprintf('%d MEL filters', nfilters)), xlabel('f')
  subplot(2,1,2)
  plot(f2, Phi(show_mel,1:nfft2));
  title(sprintf('MEL filter', show_mel)), xlabel('f')
  xlabel('f')
  pause(10);
  subplot(1,1,1)
  %if save_plots, print  2_mel_filters.pdf, end
end

plot(D(3:end))
title('Distorsion of lsee_mstftm algorithm'), xlabel('iteration')
pause(10);

if show_frame > 0   
  % Example: Plot |X|, |Xr|, etc of one frame
  nfr=100; %Show frame nfr
  infr=(nfr-1)*(interpolate_nframes+1)+1; %Corresponding interpolated

  % Plot x(t) and X(f)
  n_ini = nfr * len * shift_factor; % aprox.
  n_end = n_ini + len -1;
  t = [n_ini: n_end]/fs;
  subplot(2,1,1), plot(t, x(n_ini:n_end) .* win),       title(sprintf('x(t), frame %d', nfr)), xlabel('t') 
  subplot(2,1,2), plot(f2, 20*log10(Xm(1:nfft2, nfr))), title(sprintf('|X(f)|db, frame %d', nfr)), xlabel('f') 
  %if save_plots, print  3_x_and_X.pdf, end
  pause(10);
  subplot(1,1,1)

  % Plot X(f) and Em(f)
  plot(f2, 20*log10(Xm(1:nfft2,nfr)), fc, 20/gamma*log10(Em(:,nfr)), '*r', fc, 20/gamma*log10(Em(:,nfr)), '-r')
  title(sprintf('|X(f)|, and energy bands, frame %d', nfr)), xlabel('f') 
  %if save_plots, print  4_X_and_Eb.pdf, end
  pause(10);

  % For the moment, Em and Emq is the same: no quantization
  plot(fc, 20/gamma*log10(Em(:,nfr)),'b', fc, 20/gamma*log10(Em(:,nfr)), 'r', fc, 20/gamma*log10(Emq(:,nfr)), 'b', fc, 20/gamma*log10(Emq(:,nfr)), 'r')
  title(sprintf('Original and reconstructed energy bands, frame %d', nfr)), xlabel('f') 
  %if save_plots, print  5_Eb_and_Ebq.pdf, end
  pause(10)

  % Plot X(f) and Xmr(f)
  plot(f2, 20*log10(Xm(1:nfft2,nfr)), 'b', f2, 20*log10(Xmr(1:nfft2,infr)), 'r')
  title(sprintf('|X(f)|db and |X(f)|rec., frame', nfr)), xlabel('f') 
  %if save_plots, print  6_X_and_Xr.pdf, end
end
