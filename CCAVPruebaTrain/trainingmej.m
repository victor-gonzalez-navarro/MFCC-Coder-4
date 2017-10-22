function fmatrix = trainingmej(Wmatrixes,DCT,H,MAG,C,K,frames)
disp('Entro en la funci�n trainingmej');


%% N�mero de lapabras en el codebook (es decir, n�mero de centroides)
J=2^12;
disp('Hay un valor toal de centroides de:');
disp(J);


%% Inicializaci�n de fmatrix (F={~f1, ~f2, ... , ~fJ} que contiene por columnas los clusters

% Inicializaci�n de fmatrix destribuida unifromemente por todos los puntos Y={y1, y2,..., yK}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fmatrix = zeros(C+1,J);
for i=1:J 
    numframes = size(frames,2);
    index = floor((numframes*i)/J);
    FBE = H * MAG(1:K,index);
    fmatrix(:,i) = DCT * log( FBE );
end

% Inicializaci�n de fmatrix destribuida a las primeras ys Y={y1, y2,..., yK}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fmatrix = zeros(C+1,J);
%for i=1:J 
%    index = i;
%    FBE = H * MAG(1:K,index);
%    fmatrix(:,i) = DCT * log( FBE );
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Linea 1 del paper: inicializaci�n: n=1, T=50, et=1e8, ?=0.1,
numIter=0;
delta = 0.001; %0.1
incerror = 1e8;
T = 50;

%% Inicializaci�n de variables que se usar�n m�s alante
% Lo de K+1 es para a�adir una etiqueta para saber que y va a cada codeword
MAGaug = MAG(1:(K+1-1),:);
veczeros = zeros(1,size(MAG,2));
MAGaug = [MAGaug; veczeros];
MAGmar = zeros(K,J);
errorant = 0.0;





%% Linea 2 del paper: while n ? T,et ? ? do
while numIter<=T && incerror>=delta 
    
    %% Linia 4 del paper: reconstruir ~yj a partir de ~fj
    % Cada columna de MAGmar corresponde a ~yj=H^(-1)*exp[D^(-1)*~fj]
    disp('Linea 4 del paper');
    for j=1:J
        DCTmenos1 = inv(DCT);
        PRU = DCTmenos1 * fmatrix(:,j);
        MAGmar(:,j) = pinv(H) * exp(PRU);
    end
            
    
    %% Linia 6 del paper: asigna cada yk a un cluster Cj
    % Cj = {yk | d(yk,~yj) < d(yk,~yi), i=1,...,J; i!=j}
    disp('Linea 6 del paper');

    
    % Buscar que columna de MAG(,) [REDUCIDO, tines q ser del mismo tama�o que MAGmar] 
    % se parece mas a cualquier columna de MAGmar(,)
    numero = 1;
    for t=1:size(frames,2)
            
        if t>(30000*numero)
            disp('He pasado 50 000');
            numero = numero + 1;
        end
        
        min = 10e28;            
        v1 = MAGaug(1:K,t);
        for j=1:J
            v2 = MAGmar(1:K,j);
            if norm(v2-v1)^2 < min
                %------------------------------------------------------------------------------------------------
                W = zeros(K,K);
                for m=1:K
                     W(m,m) = Wmatrixes(m,t);
                end
                %------------------------------------------------------------------------------------------------
                
                min = (v1-v2)'* W *(v1-v2);
                MAGaug((K+1),t) = j;
            end
        end
    end
    

    %% Linia 8 del paper: actualiza el centroide ~fj :
    % Utiliza una f�rmula demostrada en el paper muy larga
    disp('Linea 8 del paper');

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % S A C A R
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %for l=1:size(frames,2)
    %    if MAGaug(K+1,l) == 17
    %        disp('PORFIN!!');
    %        l
    %    end
    %end
    %pause(20);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    for j=1:J 
        sum1 = zeros(K,K);
        sum2 = zeros(K,1);
        
        for l=1:size(frames,2)
            if MAGaug(K+1,l) == j
                %%% MEEEEEEEEEEEEJJJJJJJJJJJJJJOOOOOOOOOOOORRRRRRRRRAAAAAAARRRRRRRRR-- H es MxK y
                %%% Wmatrixes deberia ser (K-1)x(K-1)
                
                %------------------------------------------------------------------------------------------------
                W = zeros(K,K);
                for m=1:K
                     W(m,m) = Wmatrixes(m,l);
                end
                %------------------------------------------------------------------------------------------------
                
                sum1 = sum1 + W;
                sum2 = sum2 + W * MAG(1:K,l); %sum2 es un vector en vertical
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % S A C A R
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ver = (diag(sum1)==zeros(K,1));
        if ver == ones(K,1)
            for i=1:K
                sum1(i,i) = 1;
            end
            
            J = J-1;
            MAGmar(:,j)=[];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        PRUTOT = H*(inv(sum1))*sum2;
        
        fmatrix(:,j) = DCT*(log(PRUTOT)); %dim de fmatrix es Nx1 = DCT=(N*M) H(M*K) sum1(K*K) sum2(K*1)
    end
    
    
    %% Linia 10 del paper: calcula la distorsi�n total et (entre puntos y centroides):
    disp('Linea 10 del paper');

    sum3 = 0;
    for j=1:J 
        sum4 = 0;
        for l=1:size(frames,2)
            if MAGaug(K+1,l) == j
                ant1 = MAG(1:K,l);
                ant21 = inv(DCT) * fmatrix(:,j);  
                ant22 = pinv(H) * exp(ant21);
                ant = ant1 - ant22;
                
                %------------------------------------------------------------------------------------------------
                W = zeros(K,K);
                for m=1:K
                     W(m,m) = Wmatrixes(m,l);
                end
                %------------------------------------------------------------------------------------------------

                v = sqrt(W)*ant;
                res = norm(v)^2; 
                sum4 = sum4 + res;
            end
        end 
        sum3 = sum3 + sum4;
    end
    
    
    %% Linia 12 del paper: ||et^(n)-et^(n-1)||2,2
    % Calcula la variaci�n de error entre puntos y clusters de ahora y de la it. anterior
    disp('Linea 12 del paper');
    
    errornew = sum3;
    incerror = norm(errornew-errorant)^2;
    %incerror = abs(errornew-errorant);
    disp(numIter);
    disp(errornew);
    disp('J');
    disp(J);
    
    errorant = errornew;
    numIter=numIter+1;
    %numIter 
    disp('Final linea 12 del paper');
end


disp('Salgo de la funci�n trainingmej');
