close all hidden
clear all
clc

% Projet TNS TSI

%% *** Paramètres *** %%

method='butter';

% Gabarit
fg=transpose([0,0.08;0.12,0.25;0.12,0.25;0.27,0.3;0.32,0.4;0.32,0.4;0.44,0.5]);
GdBg=transpose([-55,-55;-1,-1;1,1;-20,-20;-1,-1;1,1;-55,-55]);


switch method

    case 'fir2'
%% Synthèse des filtres RIF
Npts=1024*4;
N=200;
F=2*[0,0.09,0.12,0.25,0.27,0.3,0.32,0.4,0.43,0.5];
Amp=10.^([-55,-55,0,0,-20,-20,0,0,-55,-55]);

[B,A]=fir2(N,F,Amp,Npts);
[Hrif1,Wrif1] = freqz(B,A,Npts);
Wrif1=Wrif1/(2*pi);

% Tracé
figure;
hold on
plot(Wrif1,20*log10(abs(Hrif1)))
plot(fg,GdBg,'r')
title(['Gain (dB) du filtre RIF ' method ' d' 180 'ordre ' int2str(N)]); 
xlabel('frequence normalisee'); 
ylabel('Gain (dB)'); 
legend('Filtre','Gabarit')


[phi,Wp]=phasez(B,A,Npts);

for k=1:(size(phi,1)-1)
    if phi(k+1,1)>phi(k,1)
        phi(k+1:end,1)=phi(k+1:end,1)-pi;
    end
end

TPG=-1/(2*pi)*diff(phi);

figure
plot(Wrif1,phi)
title(['Phase du filtre RIF ' method ' d' 180 'ordre ' int2str(N)]); 
xlabel('frequence normalisee'); 
ylabel('Phase (radian)'); 
legend('Filtre')


figure
plot(Wrif1(1:end-1),TPG)
t=mean(TPG);
v=var(TPG);
axis([0 0.5 t/2 3*t/2])
title(['Temps de propagation de groupe du filtre RIF ' method ' d' 180 'ordre ' int2str(N)]); 
xlabel('frequence normalisee'); 
ylabel('Temps de propagation de groupe (s)'); 
legend('Filtre') 
    
    
    
    case 'butter'
%% Synthèse des filtres RII butter
Npts=2^10;
% Filtre n°1
W11 = 0.115;
W21 = 0.255;
N1 = 14; % Ordre du premier RII
Wn1 = 2*[W11 W21];
[B1,A1] = butter(N1,Wn1);
[H1,W1] = freqz(B1,A1,Npts);
W1=W1/(2*pi);

% Filtre n°2
W12 = 0.315;
W22 = 0.405;
N2 = 10; % Ordre du deuxième RII
Wn2 = 2*[W12 W22];
[B2,A2] = butter(N2,Wn2);
[H2,W2] = freqz(B2,A2,Npts);
W2=W2/(2*pi);


    case 'cheby2'
%% Synthèse des filtres RII cheby2
Npts=2^10;
% Filtre n°1
W11 = 0.1;
W21 = 0.275;
N1 = 10; % Ordre du premier RII
Wn1 = 2*[W11 W21];
R1=55;
[B1,A1] = cheby2(N1,R1,Wn1);
[H1,W1] = freqz(B1,A1,Npts);
W1=W1/(2*pi);

% Filtre n°2
W12 = 0.30;
W22 = 0.42;
N2 = 10; % Ordre du deuxième RII
Wn2 = 2*[W12 W22];
R2=55;
[B2,A2] = cheby2(N2,R2,Wn2);
[H2,W2] = freqz(B2,A2,Npts);
W2=W2/(2*pi);


    case 'ellip'
%% Synthèse des filtres RII ellip
Npts=2^12;
% Filtre n°1
W11 = 0.11;
W21 = 0.255;
N1 = 5; % Ordre du premier RII
Wn1 = 2*[W11 W21];
Rp1=1;
Rs1=55;
[B1,A1] = ellip(N1,Rp1,Rs1,Wn1);
[H1,W1] = freqz(B1,A1,Npts);
W1=W1/(2*pi);

% Filtre n°2
W12 = 0.319;
W22 = 0.405;
N2 = 4; % Ordre du deuxième RII
Wn2 = 2*[W12 W22];
Rp2=1;
Rs2=55;
[B2,A2] = ellip(N2,Rp2,Rs2,Wn2);
[H2,W2] = freqz(B2,A2,Npts);
W2=W2/(2*pi);

end


rii=0;

if strcmp(method,'fir2')
    rii=0;
else
    rii=1;
end

rii=boolean(rii);

if rii

    % Tracé

    figure; 
    plot(W1,20*log10(abs(H1)),'b'); 
    hold on
    plot(W2,20*log10(abs(H2)),'c'); 
    plot(fg,GdBg,'r')

    axis([0 0.5 -60 2])
    title(['Synthese des RII ' method ' d' 180 'ordre ' int2str(N1) ' et ' int2str(N2) ]); 
    xlabel('f normalisee'); 
    ylabel('Gain (dB)'); 
    legend('Filtre n°1','Filtre n°2','Gabarit')
    
  
    [H1,phi1,W1p]=deroulement_phase(H1,W1);
    [H2,phi2,W2p]=deroulement_phase(H2,W2);
    
    
    TPG1=-1/(2*pi)*diff(phi1);
    TPG2=-1/(2*pi)*diff(phi2);
    
    figure
    hold on
    plot(W1,phi1,'b')
    plot(W2,phi2,'c')
    title(['Phase des RII ' method ' d' 180 'ordre ' int2str(N1) ' et ' int2str(N2) ]); 
    xlabel('frequence normalisee'); 
    ylabel('phase (radian)'); 
    legend('Filtre n°1','Filtre n°2')

    figure
    hold on
    plot(W1(1:end-1),TPG1,'b');
    plot(W2(1:end-1),TPG2,'c');
    title(['Temps de propagation de groupe des RII ' method ' d' 180 'ordre ' int2str(N1) ' et ' int2str(N2) ]); 
    xlabel('frequence normalisee'); 
    ylabel('Temps de propagation de groupe (s)'); 
    legend('Filtre n°1','Filtre n°2')

%     return
    
    % Synthèse de déphaseurs pur %
    
    
    
    t=1;
    

    ob1=2046;

    [Brif1,Arif1]=rif_dephaseur(phi1,ob1,W1p);
    [Hrif1,Wrif1] = freqz(Brif1,Arif1,Npts);
    Wrif1=Wrif1/(2*pi);
    
    
    ob2=2046;

    [Brif2,Arif2]=rif_dephaseur(phi2,ob2,W2p);
    [Hrif2,Wrif2] = freqz(Brif2,Arif2,Npts);
    Wrif2=Wrif2/(2*pi);    
    


    figure
    plot(Wrif1,20*log10(abs(Hrif1)))
    title(['Module du RIF pour ob1=' int2str(ob1)])
    
    figure
    hold on
    plot(Wrif1,(unwrap(angle(Hrif1))+phi1))
    title(['Erreur de phase pour ob1=' int2str(ob1)])
    
    
    figure
    plot(Wrif2,20*log10(abs(Hrif2)))
    title(['Module du RIF pour ob2=' int2str(ob2)])
    
    figure
    hold on
    plot(Wrif2,(unwrap(angle(Hrif2))+phi2))
    title(['Erreur de phase pour ob2=' int2str(ob2)])
    


    tracer(t,1)=ob1;
    tracer(t,2)=max(abs(unwrap(angle(Hrif1))+phi1));
    tracer(t,3)=max(20*log10(abs(Hrif1)))-min(20*log10(abs(Hrif1)));
    t=t+1;

end

figure
hold on
plot(tracer(:,1),tracer(:,2))
plot(tracer(:,1),tracer(:,3),'r')
xlabel('Ordre du filtre')
legend('Erreur de phase','Erreur de module du RIF')


