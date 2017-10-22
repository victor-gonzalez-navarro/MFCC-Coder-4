%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration and definition of mel filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename
%basefile='ona8cs';
basefile='SA000S01';


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
nfilters    = 60;
nfilters1k  = 30;

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
% C O D E R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wav_file = 'ona8.wav';  % input audio filename

%Read speech samples, sampling rate and precision from file
[ speech, fs] = audioread( wav_file );
path = fullfile('/Users/Victor/Documents/MATLAB/CCAVPruebaTrain/AudioFiles');
files = dir(path);

for fileIndex=1:length(files)
   if (~isempty(strfind(files(fileIndex).name,'wav')))
         disp(fullfile(path,files(fileIndex).name))
         [data,fs] = audioread(fullfile(path,files(fileIndex).name));
         speech = [speech ; data];
   end
end
x = speech;
%[q_mfcc Em Xm] = mfcc_coder(Phi, win, shift_factor, x, gamma);
% The coded signal is just q_mfcc. Em and Xm are just for debug and print results. 

M = size(Phi,1);
disp('Entro en ka función que pasa de señal a mfcc');
[q_mfcc, Em, Xm, DCT, Xm2] = mfcc_coder_vicnat(Phi, win, shift_factor, x, gamma,M);
disp('Esta función de mfcc_coder es nuestra y no de bonafontee');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q U A N T I Z E R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%L = 2^9;
%disp('Hay un valor toal de centroides de:');
%disp(L);
%[mVQ] = vqsplit(q_mfcc, L);
disp('Estoy usando mi quantizer del training');


K = size(Phi,2);

p = 10;
Wmatrixes = zeros(K,K);
disp('Entro en el for de Wmatrixes');
for i=1:(size(Xm,2))
    %Wmatrixes(:,:,i) = new_matrix_W(Xm(:,i), K, p);
    Wmatrixes(:,i) = new_matrix_W(Xm(:,i), K, p);
end
disp('Salgo del for de Wmatrixes');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%size(Phi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C = size(q_mfcc,1) -1;
fmatrix = trainingmej(Wmatrixes,DCT,Phi,Xm2,C,K,Xm);



