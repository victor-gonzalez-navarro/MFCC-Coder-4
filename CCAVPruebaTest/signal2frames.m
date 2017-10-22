function xm = signal2frames(xp, shift_factor, win, do_not_add_guards)
  % xm = signal2frames(x, shift_factor, win, do_not_add_guards)
  % Split and window signal x, in frames of the size = length(win).
  % The shift_factor controls the overlapping between frames: (0,1] 
  %     0 => 100% overlapp (=> infinite frames!!)
  %     1 => no overlapp
  % xm(m,n): sample 'n' of frame 'm'
  % By default, some zeros are added at the beg. and end so that signal
  % can be recovered, even if window is applied

  % ---------------------------------------------------------------------
  % Copyright (C) Antonio Bonafonte, 2014
  % Universitat Politecnica de Catalunya, Barcelona, Spain.
  %
  % Permission to copy, use, modify, sell and distribute this software
  % is granted provided this copyright notice appears in all copies.
  % This software is provided "as is" without express or implied
  % warranty, and with no claim as to its suitability for any purpose.
  % ---------------------------------------------------------------------
  len = length(win);
  shi = len * shift_factor;



  if (nargin > 3 && do_not_add_guards ~= 0)
    x = xp;
  else
    g = length(win)-shi;
    x = [zeros(g,1); xp; zeros(g,1)];
  end

  nframe = 1+ceil((length(x)-len)/shi);
  xm = zeros(len, nframe);

  for n=1:nframe
    first = round((n-1)*shi + 1);
    last = first+len-1;
    if (last <= length(x))
       wframe = x(first:last);
    else
       missing = last-length(x); 
       wframe = [x(first:end); zeros(missing,1)];
    end
    xm(:,n) = wframe .* win;
  end

end
