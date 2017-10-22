function [Phi, fc] = melfilters(N, fm, M, M1, normalize)
  % Create the matrix Phi(M x N) with M mel filters
  % fc: center frequencies
  % Usage: melfilters(N, fm, M)
  % Usage: melfilters(N, fm, M, M1)
  % Usage: melfilters(N, fm, M, M1, norm)
  % N: Number of coefficients
  % fm: sampling frequency
  % M: Number of filters
  % M1: Number of linear-spaced filters under 1k
  % normalize: 0 (def.): max. of filter 1; 1: sum coef=1  


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

  if nargin == 3
    %center frequencies
    maxmel = 1127*log(1+(fm/2)/700);
    fc=700*(exp((0:maxmel/(M-1):maxmel)/1127)-1);
  else
    fLin=1000;
    fc = 0:fLin/(M1-1):fLin;
    mel1K = 1127*log(1+1000/700);
    maxmel = 1127*log(1+(fm/2)/700);
    M2 = M - M1 + 1;
    fc2=700*(exp((mel1K:(maxmel-mel1K)/(M2-1):maxmel)/1127)-1);
    fc = [fc, fc2(2:end)];
  end
  %in samples
  fcs=1+round(fc/(fm/2) * (N/2-1));

  % reverse=1: add a reversed filter, to include F:0.5 to 1.
  % In principle, it is the same (|X| is even).

  reverse=1;

  Phi = zeros(M,N+1);
  Phi(1,fcs(1):fcs(2)) = triangle(fcs(1), fcs(1), fcs(2));
  Phi(1,:) = Phi(1,:) + reverse* fliplr(Phi(1,:));
  for i=2:M-1
    Phi(i,fcs(i-1):fcs(i+1)) = triangle(fcs(i-1), fcs(i), fcs(i+1));
    Phi(i,:) = Phi(i,:) + reverse*fliplr(Phi(i,:));
  end
  Phi(M,fcs(M-1):fcs(M)) = triangle(fcs(M-1), fcs(M), fcs(M));
  Phi(M,:) = Phi(M,:) + reverse*fliplr(Phi(M,:));

  Phi=Phi(:,1:N);

  for i=1:N
      
  end
  
  
  if (nargin >= 5 && normalize ~= 0)
      vec = sum(Phi,2);
      for i=1:M
        Phi(i,:) = Phi(i,:) / vec(i); %?????????????????????????????????????????????????????????????????????
      end
  end
end




function y = triangle(a, b, c)
  y = zeros(c-a+1,1);
  if (b > a)
    i = a:b;
    y(i-a+1) = (i-a)/(b-a);
  end
  if (c > b)
    i = b:c;
    y(i-a+1) = (c-i)/(c-b);
  end
end
