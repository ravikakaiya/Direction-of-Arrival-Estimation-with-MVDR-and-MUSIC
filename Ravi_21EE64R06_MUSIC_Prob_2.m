%----------------------------------
%Author: Ravi Kakaiya
%Roll No : 21EE64R06
%----------------------------------
%RTSP LAB PROBLEM : 2
%--------------------------------------------------------------------------
%Part 1 :Simulation of the given System and generating data matrix 
%--------------------------------------------------------------------------
%Given Parameters:
%--------------------------------------------------------------------------
%Frequency
f=(2*(10^9));
%wavelength 
lamda=(3*(10^8)/f);
d=lamda/2;T=4500;
%number of Nodes
Nr=10;e=exp(1);
%angls of signals in rad.
DOA=[-10,20];
K=length(DOA);
%-------------------------------------------------
% Steering vector matrix.
A=exp(-1i*2*pi*(d/lamda)*(0:Nr-1)'*sind(DOA));
%Generating random signal
sig=randn(K,T); 
sig=e.^(1i*2*pi*f*sig);
%-------------------------------------------------
%Ploting generated Signal 
figure;
subplot 211
plot(1:length(sig(1,:)),abs(sig));
title('Random Signal Magnitude plot')
ylabel('Magnitude')
xlabel('Number of Datasnapshot')
p1= angle(sig);
subplot 212
plot(1:length(sig(1,:)),p1);
ylabel('Magnitude')
xlabel('Number of Datasnapshot')
title('Random Signal phase plot')
%-------------------------------------------------
%Uncorrelated white noise
n=(randn(Nr,T)); 
%Generating data matrix
X=A*sig+n; 
%-------------------------------------------------
%Observing White Noise and Generated Data Matrix
px = abs(X(:,1:10))
figure;
subplot 121
imagesc(px)
title('Data matrix (Fraction Part)')
subplot 122
imagesc(abs(n(:,1:10)))
title('White Noise (Fraction Part)')
%--------------------------------------------------------------------------
%Part 2 :Estimating Peak for Direction Of Arrivals
%--------------------------------------------------------------------------
%Calculating Covariance matrix
Rx=X*ctranspose(X); 
%Plotting Covariance Matrix
figure;
imagesc(abs(Rx))
xlabel('Number of sensor Nodes')
ylabel('Number of sensor Nodes')
title("Visualizing Covariance Matrix")
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
subplot 211
plot(theta,abs(music_spect))
title('MUSIC Spectrum')
xlabel('Angle in degrees')
ylabel('Magnitude')
%--------------------------------------------------------------------------
%Plotting Bars for Peak Frequencies: 
ang = zeros(length(theta),1);
for i=1:length(theta)
    for j = 1:length(locs)
    if find(i==locas(locs(j)))
        ang(i)=pks(j);
    end
    end
end
subplot 212
bar(theta,ang,10,'r')
title('Peaks For Estimated theta')
xlabel('Angle in degrees')
ylabel('Magnitude')
%--------------------------------------------------------------------------