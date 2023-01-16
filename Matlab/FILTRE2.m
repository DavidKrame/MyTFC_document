% %R�cup�ration du son

%D�finition de la p�riode d'�chantillonnage
% fs=44100;
% %D�finition du nombre de voies d'enregistrement
% noc=1;
% %D�finition du nombre de bits
% nob=16;
% %Enregistrement du signal audio
% recObj=audiorecorder(fs,nob,noc);
% record(recObj);
% pause(40);
% stop(recObj);
% %Lancement du signal enregistr�
% play(recObj);
% %Formatage du signal enregistr� sous forme vectoriel
% myrecordings=getaudiodata(recObj);
% %Repr�sentation du signal enregistr�
% plot(myrecordings);
% %Fin de la partie prise de son

%D�but de l'algorithme et lancement du chronom�trage
tic;
%Initialisations et d�claration des constantes

w=zeros(1,64);
zero=zeros(1,64);
%Compl�ment de W par des zeros pour atteindre la longueur r�quise
wk=[w,zero]';
%D�claration des constantes "beta", "delta", "rho" et "mu".
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

%D�but du traitement
for i=1:20850
    %S�lection particulier du premier bloc d'�chantillons
    if i==1
        sign=(myrecordings(1:128))';
        sign2=(myrecordings(128:-1:1));
    %S�lection g�n�ral de tous les blocs    
    elseif i>=2
        sign=(myrecordings(65+64*(i-2):192+64*(i-2)))';
        sign2=(myrecordings((192+64*(i-2)):-1:(65+64*(i-2))));
    end
    
    %Transformation du signal d'entr�e et du vecteur W
    SIG=fft(sign);
    W=(fft(wk))';
    %Calcul du produit terme � terme des deux vecteurs
    Y=SIG.*W;
    %Transformation inverse pour retrouver le produit de convolution
    YW=ifft(Y);
    %S�lection des termes utiles recherch�s
    yw=YW(65:128);
    %Calcul du vecteur d'erreurs
    err=(myrecordings(1+64*(i-1):64*i))'-(yw);
    %Compl�ment du vecteur par des zeros puis transformation
    ERR=[err,zero];
    errF=fft(ERR);
    %Transformation du vecteur image
    SIG2=fft(sign2);
    %R�alisation de l'�tape.... de l'algorithme
    ERsign=(errF)'.*(SIG2);
    erSIGN=ifft(ERsign);
    es=erSIGN(65:128);
    %R�alisation de l'�tape ..... (Calcul du vecteur d'autocorr�lation)
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
    %R�alisation de l'�tape....
    GkRss=Gk*(Rss);
    %R�alisation de l'�tape....
    GkEsi2=Gk*es;
    
    %R�initialisation des coefficients du filtre....
    
    for j=1:64
        %Algorithme BPNLMS
        w1(j)=(w(j)+((mu*GkEsi2(j))/(GkRss(j)+b)));
        %Algorithme BNLMS
        w2(j)=(w(j)+((mu*es(j))/(Rss(j)+b)));
        
        %R�initialisation des coefficients du filtre par BPNLMS++
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
%Fin du chronom�trage
toc;