function [W] = new_matrix_W(frame, K, p)

%W = zeros(K,K);
W = zeros(K,1);
aprim = lpc(frame,p);
a = -aprim(2:p+1);


%for m=1:(K-1)
for m=1:K
    var_num = 0;
    var_den = 0;
    w = 2*pi*m/(K*2-2);
    z = exp(1i*w);
    for n=1:p
        var_num = var_num + a(n)*((0.9)^n)*z^(-n);
        var_den = var_den + a(n)*((0.5)^n)*z^(-n);
    end
    %W(m,m) = abs((1 - var_num)/(1 - var_den));
    W(m,1) = abs((1 - var_num)/(1 - var_den));
end

