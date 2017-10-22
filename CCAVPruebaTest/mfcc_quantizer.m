function [fin] = mfcc_quantizer(q_mfcc,mVQ)

 %Linia 6 del paper: asigna cada yk a un cluster Cj
    % Cj = {yk | d(yk,~yj) < d(yk,~yi), i=1,...,J; i!=j}
    
    % Buscar que columna de MAG(,) [REDUCIDO, tines q ser del mismo tamaño que MAGmar] 
    % se parece mas a cualquier columna de MAGmar(,)
    xtest = q_mfcc;
    K = (size(xtest,2));
    J = (size(mVQ,2)); 
    
    xtestAug = [xtest;zeros(1,K)];
    
    num = (size(xtest,1)); 
    
    for t=1:K
        min = 10e28;            
        v1 = xtestAug(1:num,t);
        for j=1:J
            v2 = mVQ(1:num,j);
            if norm(v2-v1)^2 < min
                min = norm(v2-v1)^2;
                xtestAug((num+1),t) = j;
            end
        end
    end
    
    for t=1:K
        var = xtestAug(num+1,t);
        xtest(1:num,t) = mVQ(:,var);
    end
    
 fin = xtest;

    
    
