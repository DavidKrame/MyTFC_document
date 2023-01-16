% %Récupération du son

%Définition de la période d'échantillonnage
% fs=44100;
% %Définition du nombre de voies d'enregistrement
% noc=1;
% %Définition du nombre de bits
% nob=16;
% %Enregistrement du signal audio
% recObj=audiorecorder(fs,nob,noc);
% record(recObj);
% pause(40);
% stop(recObj);
% %Lancement du signal enregistré
% play(recObj);
% %Formatage du signal enregistré sous forme vectoriel
% myrecordings=getaudiodata(recObj);
% %Représentation du signal enregistré
% plot(myrecordings);
% %Fin de la partie prise de son

%Début de l'algorithme et lancement du chronométrage
tic;
%Initialisations et déclaration des constantes

w=zeros(1,64);
zero=zeros(1,64);
%Complément de W par des zeros pour atteindre la longueur réquise
wk=[w,zero]';
%Déclaration des constantes "beta", "delta", "rho" et "mu".
b=100;
D=(1/64);
rho=(10^(-2));
mu=(0.8);
%Initialisation des W pour BNLMS puis BPNLMS.
w1=(ones(1,64))';
w2=(ones(1,64))';
%Initialisation des vecteurs utiles au calcul de Gk.
gam=zeros(1,64);
gk=zeros(1,64);

%Début du traitement
for i=1:20850
    %Sélection particulier du premier bloc d'échantillons
    if i==1
        sign=(myrecordings(1:128))';
        sign2=(myrecordings(128:-1:1));
    %Sélection général de tous les blocs    
    elseif i>=2
        sign=(myrecordings(65+64*(i-2):192+64*(i-2)))';
        sign2=(myrecordings((192+64*(i-2)):-1:(65+64*(i-2))));
    end
    
    %Transformation du signal d'entrée et du vecteur W
    SIG=fft(sign);
    W=(fft(wk))';
    %Calcul du produit terme à terme des deux vecteurs
    Y=SIG.*W;
    %Transformation inverse pour retrouver le produit de convolution
    YW=ifft(Y);
    %Sélection des termes utiles recherchés
    yw=YW(65:128);
    %Calcul du vecteur d'erreurs
    err=(myrecordings(1+64*(i-1):64*i))'-(yw);
    %Complément du vecteur par des zeros puis transformation
    ERR=[err,zero];
    errF=fft(ERR);
    %Transformation du vecteur image
    SIG2=fft(sign2);
    %Réalisation de l'étape.... de l'algorithme
    ERsign=(errF)'.*(SIG2);
    erSIGN=ifft(ERsign);
    es=erSIGN(65:128);
    %Réalisation de l'étape ..... (Calcul du vecteur d'autocorrélation)
    RSS=SIG2.*(SIG)';
    RSs=ifft(RSS);
    Rss=RSs(65:128);
    
    %Calcul de .... pour le vecteur Gk
    Vk=max([D,abs(w)]);
    %Calcul du vecteur Gk
    for k=1:64
        gam(k)=max(((rho)*Vk),(abs(w(k))));
    end
    for M=1:64
        gk(M)=((gam(M))/(mean(gam)));
    end
    Gk=diag(gk);
    %Réalisation de l'étape....
    GkRss=Gk*(Rss);
    %Réalisation de l'étape....
    GkEsi2=Gk*es;
    
    %Réinitialisation des coefficients du filtre....
    
    for j=1:64
        %Algorithme BPNLMS
        w1(j)=(w(j)+((mu*GkEsi2(j))/(GkRss(j)+b)));
        %Algorithme BNLMS
        w2(j)=(w(j)+((mu*es(j))/(Rss(j)+b)));
        
        %Réinitialisation des coefficients du filtre par BPNLMS++
        %Coefficients impairs par BPNLMS
        if mod(j,2)==1
            w(j)=w1(j);
        %Coefficients pairs par BNLMS    
        elseif mod(j,2)==0
            w(j)=w2(j);
        end
    end
end
% w=((ones(1,64))*(13.5)-(w))/(13.5);
%Fin du chronométrage
toc;