% Projet TNS TSI

%% Synthèse des filtres RII butter
clc; clear all; close all;
Npts=2^10;
% Filtre n°1
W11 = 0.115;
W21 = 0.255;
N1 = 14; % Ordre du premier RII
Wn1 = 2*[W11 W21];
[B1,A1] = butter(N1,Wn1);
[H1,W1] = freqz(B1,A1,Npts);
[phase1,Wphase1] = phasez(B1,A1,Npts);
W1=W1/(2*pi);
Wphase1=Wphase1/(2*pi);
phase1=phase1*180/pi;

% Filtre n°2
W12 = 0.315;
W22 = 0.405;
N2 = 10; % Ordre du deuxième RII
Wn2 = 2*[W12 W22];
[B2,A2] = butter(N2,Wn2);
[H2,W2] = freqz(B2,A2,Npts);
[phase2,Wphase2] = phasez(B2,A2,Npts);
W2=W2/(2*pi);
Wphase2=Wphase2/(2*pi);
phase2=phase2*180/pi;

% Gabarit
fgab=[0,0.08;0.12,0.25;0.12,0.25;0.27,0.3;0.32,0.4;0.32,0.4;0.44,0.5]';
GdBgab=[-55,-55;1,1;-1,-1;-20,-20;1,1;-1,-1;-55,-55]';

% Tracé
figure; % Synthese des RII
plot(W1,20*log10(abs(H1)),'b'); 
hold on
plot(W2,20*log10(abs(H2)),'b'); 
plot(fgab,GdBgab,'r')
axis([0 0.5 -60 2])
title(['Synthese des RII d ordre ' int2str(N1) ' et ' int2str(N2) ]); 
xlabel('f normalisee'); 
ylabel('Gain');
figure; % Phase des filtres
subplot(2,1,1);
plot(Wphase1,phase1,'b');
axis([0 0.5 min(phase1) max(phase1)]);
subplot(2,1,2);
plot(Wphase2,phase2,'b');
axis([0 0.5 min(phase2) max(phase2)]);
xlabel('f normalisee'); 
ylabel('Phase');

%% Synthèse des filtres RII cheby2
clc; clear all; close all;
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
W12 = 0.29;
W22 = 0.42;
N2 = 8; % Ordre du deuxième RII
Wn2 = 2*[W12 W22];
R2=55;
[B2,A2] = cheby2(N2,R2,Wn2);
[H2,W2] = freqz(B2,A2,Npts);
W2=W2/(2*pi);

% Gabarit
fgab=[0,0.08;0.12,0.25;0.12,0.25;0.27,0.3;0.32,0.4;0.32,0.4;0.44,0.5]';
GdBgab=[-55,-55;1,1;-1,-1;-20,-20;1,1;-1,-1;-55,-55]';

% Tracé
figure; 
plot(W1,20*log10(abs(H1)),'b'); 
hold on
plot(W2,20*log10(abs(H2)),'b'); 
plot(fgab,GdBgab,'r')
axis([0 0.5 -60 2])
title(['Synthese des RII d ordre ' int2str(N1) ' et ' int2str(N2) ]); 
xlabel('f normalisee'); 
ylabel('Gain'); 
