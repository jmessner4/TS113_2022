% % Messner Julie & Penot Victorine

clear ;     % Efface les variables de l ’ environnement de travail
close all ; % Ferme les figures ouvertes
clc ;       % Efface la console

% % Initialisation des paramètres

fe = 1e4 ;  % Fréquence d ’ échantillonnage
M = 4;      % Nombre de symboles dans la modulation
n_b = log2 ( M ) ;  % Nombre de bits par symboles

% ... autres paramètres

bits = M*n_b;
Ds = 1e3;
Ts = 1e-3;
Fse = Ts*fe;
g=@(t) t>0 & t<Ts;
t = linspace(-1e-2,1e-2,1000);
gt = ones(Fse,1);



Ns = 1000;
Nfft = 512;

% Energie du filtre de mise en forme
Eg =  1*Fse;

sigA2 = 0;                                % Variance théorique des symboles
eb_n0_dB = 0:0.5:10;                        % Liste des Eb / N0 en dB
eb_n0 = 10.^( eb_n0_dB /10) ;               % Liste des Eb / N0
sigma2 = sigA2 * Eg ./ ( n_b * eb_n0 ) ;    % Variance du bruit complexe en bande de base
TEB = zeros ( size ( eb_n0 ) ) ;            % Tableau des TEB ( résultats )
Pb = qfunc ( sqrt (2* eb_n0 ) ) ;           % Tableau des probabilités d ’erreurs théoriques
ga=gt;                                      %Filtre adapté

%Symboles
c1= 1/sqrt(2) + 1j/sqrt(2);
c2= -1/sqrt(2) + 1j/sqrt(2);
c3= -1/sqrt(2) - 1j/sqrt(2);
c4= 1/sqrt(2) - 1j/sqrt(2);


for i = 1: length ( eb_n0 )
    error_cnt = 0;
    bit_cnt = 0;
    while error_cnt < 100
        % % Émetteur
        X = randi([0,1],10000,1);
        %Association bits->symbole
        Ss = zeros(5000,1);
        k=1;
        for n=1:2:length(X) 
            if X(n) == 0 && X(n+1) == 0
                Ss(k)= c1;
                k = k+1;
            elseif X(n) == 0 && X(n+1) == 1
                Ss(k)= c2;
                k = k+1;
            elseif X(n) == 1 && X(n+1) == 1
                Ss(k)= c3;
                k = k+1;
            elseif X(n) == 1 && X(n+1) == 0
                Ss(k)= c4;
                k = k+1;
            end
        end
                Ssech = upsample(Ss,Fse);


        % Filtre de mise en forme

        Sl = conv(Ssech,gt);
        Sllen=length(Sl);
        % % Canal
        nl = sqrt (sigma2 (i)/2)*(randn(size(Sl))+1j*randn(size(Sl))) ; % Génération du bruit blanc gaussien complexe
        yl=Sl+nl;
        % % Récepteur
        
        Rl = conv(yl,ga);
        Rllen=length(Rl);

        %Echantillonnage

        Rl2 = zeros(5000,1);
        k=1;
        for p=Fse:Fse:length(Rl)-Fse
            Rl2(k)=Rl(p);
            k =k+1;
        end
        Rl2len= length(Rl2);
        
        % % Décodeur en ligne
        %Décision 

        An=zeros(Rl2len,1);

        for n=1:Rl2len
            re=real(Rl2(n));
            im=imag(Rl2(n));
            if im>=0
                if re>=0
                    An(n)=0;
                elseif re<0
                    An(n)=1;
                end
            elseif im<0
                if re>=0
                    An(n)=2;
                elseif re<0
                    An(i)=3;
                end
            end
        end

        Anlen=length(An);

        %Association Symbole bit 

        D=zeros(10000,1);

        for n=1:Anlen 
            if An(n)==0
                D(2*n-1)=0;
                D(2*n)=0;
            elseif An(i)==1
                D(2*n-1)= 0;
                D(2*n)=1;
            elseif An(n)==3
                D(2*n-1)=1;
                D(2*n)=1;
            elseif An(n)==2
                D(2*n-1)=1;
                D(2*n)=0;
            end
        end
        %error_cnt = error_cnt + sum(D~=X);
        %bit_cnt = bit_cnt + length(D);
        for n=1:length(D)
             if D(n)~=X(n)
                 error_cnt = error_cnt + 1; % incrémenter le compteur d ’ erreurs
            end
            %bit_cnt = bit_cnt + 1; % incrémenter le compteurde bits envoyés 
        end
        TEB ( i ) = error_cnt / 10000 ; 
    end
end


% % Affichage des résultats

figure,

% % Mapping de gray
semilogy(eb_n0_dB,Pb);

semilogy(eb_n0_dB,TEB),
title('Taux de l erreur binaire en fonction de eb_n0_dB')
