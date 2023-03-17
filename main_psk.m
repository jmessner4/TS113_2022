%% Messner Julie & Penot Victorine
clear;      % Efface  les  variables  de l’environnement  detravail
close  all; % Ferme  les  figures  ouvertes
clc;         % Efface  la  console

%% Initialisation  des  paramètres
fe = 1e4; % Fréquence d’échantillonnage
M = 4; % Nombre  de  symboles  dans la  modulation
n_b = log2(M); % Nombre  de bits  par  symboles
% ...  autres  paramètres
bits = M*n_b;
Ds = 1e3;
Ts = 1e-3;
Fs=1/Ts;
Fse = Ts*fe;
%g=@(t) t>0 & t<Ts;
%t = linspace(-1e-2,1e-2,1000);
gt = ones(Fse,1);

X = randi([0,1],1000,1);

Ns = 1000;
Nfft = 512;

%% Émetteur

%Association bits->symbole
%Symboles
c1= 1/sqrt(2) + 1j/sqrt(2);
c2= -1/sqrt(2) + 1j/sqrt(2);
c3= -1/sqrt(2) - 1j/sqrt(2);
c4= 1/sqrt(2) - 1j/sqrt(2);
Ss = zeros(500,1);
k=1;
for i=1:2:1000 
    if X(i) == 0 && X(i+1) == 0
        Ss(k)= c1;
        k = k+1;
    elseif X(i) == 0 && X(i+1) == 1
        Ss(k)= c2;
        k = k+1;
    elseif X(i) == 1 && X(i+1) == 1
        Ss(k)= c3;
        k = k+1;
    elseif X(i) == 1 && X(i+1) == 0
        Ss(k)= c4;
        k = k+1;
    end
end

Ssech = upsample(Ss,Fse);


% Filtre de mise en forme

Sl = conv(Ssech,gt);
Sllen=length(Sl);


%% Récepteur
%Symboles
c1= 1/sqrt(2) + 1j/sqrt(2);
c2= -1/sqrt(2) + 1j/sqrt(2);
c3= -1/sqrt(2) - 1j/sqrt(2);
c4= 1/sqrt(2) - 1j/sqrt(2);
Tabsymboles=[c1, c2, c3, c4];


%Filtre de réception

% for i=Fse:2*Fse-1
%     ga(i,1) = 1;
% end

ga=gt;

Rl = conv(Sl,ga);
Rllen=length(Rl);

%Echantillonnage

Rl2 = zeros(500,1);
k=1;
for i=Fse:Fse:length(Rl)-Fse
    Rl2(k)=Rl(i);
    k =k+1;
end
Rl2len= length(Rl2);



%Décision

An=zeros(Rl2len,1);

for i=1:Rl2len
    re=real(Rl2(i));
    im=imag(Rl2(i));
    if im>=0
        if re>=0
            An(i)=0;
        elseif re<0
            An(i)=1;
        end
    elseif im<0
        if re>=0
            An(i)=2;
        elseif re<0
            An(i)=3;
        end
    end
end

Anlen=length(An);

%Association bit Symbole

D=zeros(1000,1);

for i=1:Anlen 
    if An(i)==0
        D(2*i-1)=0;
        D(2*i)=0;
    elseif An(i)==1
        D(2*i-1)= 0;
        D(2*i)=1;
    elseif An(i)==3
        D(2*i-1)=1;
        D(2*i)=1;
    elseif An(i)==2
        D(2*i-1)=1;
        D(2*i)=0;
    end
end

erreur = 0;
for i=1:length(D)
    if X(i) ~= D(i)
        erreur = erreur+1;
    end
end

erreur = erreur/length(D);

[pxx, freq] =pwelch(Sl,ones(1,Nfft),0,Nfft,fe,'centered');

%Calcul théorique de la DSP

Gf=fft(gt);
G = sinc(freq*Ts).*exp(-i*pi*freq*Ts);
R0=1     ;   %Seul moment où l'espèrance n'est pas nulle
DSPth=((abs(Gf).*abs(Gf))/Ts)*R0*exp(0);
DSPth2=Ts * abs(sinc(Ts*freq).^2);

%% Affichage  des  résultats

figure,
hold on;
axes()
plot(Sl,'*');
title('Représentation des symboles émis');
xlabel('Réels')
ylabel('Imaginaires')


%Rlred=Rl(1:10*Ts-1/fe);
figure,
hold on;
axes()
plot(Rl,'*');
title('Représentation des symboles reçus');
xlabel('Réels')
ylabel('Imaginaires')

%Periodogramme

figure,
semilogy(freq,pxx);
title('Périodogramme');

%DSP théorique
hold on,
semilogy(freq,DSPth2);
title('DSP théorique');
