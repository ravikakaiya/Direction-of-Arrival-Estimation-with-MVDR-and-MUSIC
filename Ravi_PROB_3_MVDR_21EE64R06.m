%----------------------------------
%Author: Ravi Kakaiya
%Roll No : 21EE64R06
%----------------------------------
%RTSP LAB PROBLEM : 3
%--------------------------------------------------------------------------
%Part 1 :Simulation of the given System and generating data matrix 
%--------------------------------------------------------------------------
%Given Parameters:
%--------------------------------------------------------------------------   
clear all
clc
close all
DOA = [-10 20];
%DOA = [-50 -30 -10 0 10 30 50];      %Direction of arrival (Degree)
T   = 4500;         %Snapshots (or Samples)
K   = length(DOA); %The number of signal source(or traget) 
Nr  = 10;          %Number of receiver's antennas 
lamda = 150;      %Wavelength 
d   = lamda/2;%Receiver's antennas spacing
%d = lamda/5 ; 
SNR = 1;           %Signal to Noise Ratio (dB)
%--------------------------------------------------------------------------
% Steering vector matrix.
A = zeros(Nr,K); 
for k=1:K 
    A(:,k) = exp(-1j*2*pi*d*sind(DOA(k))*(0:Nr-1)'/lamda); 
end 
Vj = diag(sqrt((10.^(SNR/10))/2));
s = Vj* ( randn(K,T) + 1j*randn(K,T) );
noise = sqrt(1/2)*(randn(Nr,T)+1j*randn(Nr,T));
%--------------------------------------------------------------------------
%Ploting generated Signal 
figure;
subplot 211
plot(1:length(s(1,:)),abs(s));
title('Random Signal Magnitude plot')
ylabel('Magnitude')
xlabel('Number of Datasnapshot')
p1= angle(s);
subplot 212
plot(1:length(s(1,:)),p1);
ylabel('Magnitude')
xlabel('Number of Datasnapshot')
title('Random Signal phase plot')
X = A*s; 
X = X+noise;      %Insert Additive White Gaussain Noise (AWGN) 
%--------------------------------------------------------------------------
%Observing White Noise and Generated Data Matrix
px = abs(X(:,1:10))
figure;
subplot 121
imagesc(px)
title('Data matrix (Fraction Part)')
subplot 122
imagesc(abs(noise(:,1:10)))
title('White Noise (Fraction Part)')
%--------------------------------------------------------------------------
%Part 2 :Estimating Peak for Direction Of Arrivals
%--------------------------------------------------------------------------
%% MVDR (Capon)
%--------------------------------------------------------------------------
theta = -90:0.05:90;       %Grid points of Peak Search 
Rx = cov(X');                     %Data covarivance matrix 
Cov_Inv = Rx^(-1); %Inverse of covariance matrix
%Plotting Covariance Matrix
figure;
subplot 121
imagesc(abs(Rx))
xlabel('Number of sensor Nodes')
ylabel('Number of sensor Nodes')
title("Visualizing Covariance Matrix")
subplot 122
imagesc(abs(Cov_Inv))
xlabel('Number of sensor Nodes')
ylabel('Number of sensor Nodes')
title("Visualizing Inverse Covariance Matrix")
%--------------------------------------------------------------------------
for i=1:length(theta) 
    SS = zeros(Nr,1); 
    SS = exp(-1j*2*pi*d*(0:Nr-1)'*sind(theta(i))/lamda);
    %MVDR Spectrum
    PP = SS'*Cov_Inv*SS;
    P_MVDR(i) = 1/ PP; 
end
P_MVDR = real(10*log10(P_MVDR)); %Spatial Spectrum function
[pks,locs] = findpeaks(P_MVDR,theta,'SortStr','descend','Annotate','extents');
MVDR_Estim = sort(locs(1:K))
%--------------------------------------------------------------------------
%Plotting Spatial Spectrum 
figure
plot(theta,P_MVDR,'-b',locs(1:K),pks(1:K),'r*'); hold on
xlim([min(theta) max(theta)])

xlabel('Angle \theta (degree)'); ylabel('Spatial Power Spectrum P(\theta) (dB)') 
title('DOA estimation based on MVDR algorithm ') 
xlim([min(theta) max(theta)])

grid on
%--------------------------------------------------------------------------
%% Conventional Beamformer 

for i=1:length(theta) 
    SS = zeros(Nr,1); 
    SS = exp(-1j*2*pi*d*(0:Nr-1)'*sind(theta(i))/lamda);
    %Beamformer Spectrum
    PP = SS'*Rx*SS;
    PP1 = SS'*SS;
    P_Conv_Beam(i) = PP/PP1; 
end
P_Conv_Beam = real(10*log10(P_Conv_Beam)); %Spatial Spectrum function
[pks,locs] = findpeaks(P_Conv_Beam,theta,'SortStr','descend','Annotate','extents');
Conv_beam_Estim = sort(locs(1:K))
%--------------------------------------------------------------------------
%Plotting Spatial Spectrum 
figure
plot(theta,P_Conv_Beam,'-b',locs(1:K),pks(1:K),'r*'); hold on
xlim([min(theta) max(theta)])

xlabel('Angle \theta (degree)'); ylabel('Spatial Power Spectrum P(\theta) (dB)') 
title('DOA estimation based on Conventional Beamformer algorithm ') 
xlim([min(theta) max(theta)])

grid on

%% MUSIC
%Compute eigendecomposition of covariance matrix
[E ,D]=eig(Rx); 
%Find p largest eigenvalues corresponding to signal space
[D,I]=sort(diag(D),1,'descend'); 
%Sorting the eigenvectors 
E=E (:,I);
%Eigenvectors corresponding to Signal-space
Es=E (:,1:K);
%Eigenvectors corresponding to nullspace(noise)  
En=E(:,K+1:Nr); 
%--------------------------------------------------------------------------
% Define theta for MUSIC spectrum
theta=(-90:0.1:90);
music_spect=zeros(length(theta),1);
%Compute steering vectors corresponding values in theta
a1=exp(-1i*2*pi*(d/lamda)*(0:Nr-1)'*sind(theta));
%Computing MUSIC spectrum
for k=1:length(theta)
music_spect(k)=1/(ctranspose(a1(:,k))*(En*ctranspose(En))*a1(:,k));
end
%--------------------------------------------------------------------------
%Finding (p) no. of peaks in MUSIC spectrum 
[vals locas] = findpeaks(abs(music_spect));
[pks locs] = maxk(vals,K);
%Plotting MUSIC spectrum 
figure;
plot(theta,abs(music_spect))
title('MUSIC Spectrum')
xlabel('Angle in degrees')
ylabel('Magnitude')
%--------------------------------------------------------------------------