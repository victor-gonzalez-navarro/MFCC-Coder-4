function x = frames2signal(xm, shift_factor, winSyn, do_not_remove_guards, winNorm)
  % x = frames2signal(xm, shift_factor, winSyn, do_not_remove_guards, winNorm)
  % Recontruct the signal from the overlapped frames.
  % If winSyn is provided, frames xm(m,:) are windowed before the sum.
  % If winNorm is provided, signal is normalized using the sum of overlapped winN
  % The shift_factor controls the overlapping between frames: (0,1] 
  %     0 => 100% overlapp (=> inf. frames!!)
  %     1 => no overlapp
  % xm(m,n): sample 'n' of frame 'm'
  %
  % => x[N*m+n] = sum_m{xm(m,n)}
  % 
  % By default, the beg./end guards added by signal2frames are deleted. 
  % But if do_not_remove_guards is !=0, then no sample is removed


  % ---------------------------------------------------------------------
  % Copyright (C) Antonio Bonafonte, 2014
  % Universitat Politecnica de Catalunya, Barcelona, Spain.
  %
  % Permission to copy, use, modify, sell and distribute this software
  % is granted provided this copyright notice appears in all copies.
  % This software is provided "as is" without express or implied
  % warranty, and with no claim as to its suitability for any purpose.
  % ---------------------------------------------------------------------

  len = size(xm,1);
  nframes = size(xm,2);
  shi = len * shift_factor;
  total_len = round((nframes-1)*shi + len);
  x = zeros(total_len,1);
  if (nargin > 4) norm = zeros(size(x)); end
  for n=1:nframes
    first = round((n-1)*shi + 1);
    wframe = xm(:,n);
    if (nargin > 2) wframe = wframe .* winSyn; end
    x(first:first+len-1) = x(first:first+len-1) + wframe;
    if (nargin > 4) 
      norm(first:first+len-1) = norm(first:first+len-1) + winNorm;
    end
  end

  if (nargin > 4) 
    x = x ./ norm;
  end

  if (nargin <= 3 || do_not_remove_guards == 0)
    g = len-shi;
    x = x(g+1:end-g);
  end
end
