% %Récupération du son
% 
% %Définition de la période d'échantillonnage
% fs=44100;
% %Définition du nombre de voies d'enregistrement
% noc=1;
% %Définition du nombre de bits
% nob=64;
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

yk=ones(1,64);
%Déclaration des constantes "beta", "delta", "rho" et "mu".
b=100;
D=(1/64);
rho=0.01;
mu=0.8;
%Initialisation des W pour BNLMS puis BPNLMS.
w1=(ones(1,64))';
w2=(ones(1,64))';
%Initialisation des vecteurs utiles au calcul de Gk.
gam=zeros(1,64);
gk=zeros(1,64);
%Initialisation des vecteurs d'erreur et de W.
err=(1:64000);
w=zeros(1,64);

%Début du traitement
for i=1:20850
    %Récupération du vecteur d'entrée à chaque étape i
    sign=(myrecordings((64+(i-1)):-1:(1+(i-1))))';
    
    %Calcul de .... pour le vecteur Gk
    Vk=max([D,abs(w)]);
    %Calcul du vecteur Gk
    for k=1:64
        gam(k)=max(((rho)*Vk),(abs(w(k))));
    end
    for M=1:64
        gk(M)=((gam(M))/(mean(gam)));
    end
    Gk=1*diag(gk);
    %Calcul des éléments du vecteur d'erreur en gérant le début
    for s=1+(i-1):64+(i-1)
        sign2=(myrecordings((64+(s-1)):-1:(1+(s-1))))';
        yk(s) = (sign2)*w';
        err(s+63)=((myrecordings(64+(s-1))-yk(s)));
    end
    
    %Début du processus de réinitialisation du filtre
    
    %Algorithme PNLMS (Calcul du vecteur W1)
        w1=(w'+((mu*Gk*((sign)')*err(i+63))/(((sign)*Gk*((sign)')+b))));
    %Algorithme NLMS (Calcul du vecteur W2)    
        w2=(w'+((mu*((sign)')*err(i+63))/(((sign)*((sign)')+b))));
        
        
    %Réinitialisation des coefficients du filtre par PNLMS++
    
    for k=1:64
        %Coefficients impairs par PNLMS
        if mod(k,2)==1
          w(k)=w1(k);
        %Coefficients pairs par NLMS
        elseif mod(k,2)==0
            w(k)=w2(k);
        end
    end
end
%Fin du chronométrage
toc;