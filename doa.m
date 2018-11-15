clear all;
%---------------------SETUP MODEL------------------------------------------
%uniform linear array of 11 sensors with isotropic power kept at 0.5 lambda

c = physconst('LightSpeed');
fc = physconst('LightSpeed');              % Operating frequency
lambda = c/fc;  %wavelength of 1 m
uniform_distance = lambda/2;

ula = phased.ULA('NumElements',11,'ElementSpacing',uniform_distance);

%---------------------signal model-----------------------------------------
ang1 = [-60; 0];          % First signal
ang2 = [-47; 0];          % Second signal
ang3 = [-34; 0];
ang4 = [-20; 0];
ang5 = [-7; 0];
ang6 = [6; 0];
ang7 = [20; 0];
ang8 = [33; 0];
ang9 = [46; 0];
ang10 = [60; 0];

%signal model with 10 sources at equally spaced between -60 and +60
angs = [ang1 ang2 ang3 ang4 ang5 ang6 ang7 ang8 ang9 ang10 ];
ang = [ang1(1,1) ang2(1,1) ang3(1,1) ang4(1,1) ang5(1,1) ang6(1,1)...
    ang7(1,1) ang8(1,1) ang9(1,1) ang10(1,1)]/180*pi;

% Another set of angles for another experiment
ang1 = [-60; 0];          % First signal
ang2 = [-50; 0];          % Second signal
ang3 = [-34; 0];
ang4 = [-31; 0];
ang5 = [-20; 0];
ang6 = [-5; 0];
ang7 = [-8; 0];
ang8 = [5; 0];
ang9 = [10; 0];
ang10 = [25; 0];
ang11 = [41; 0];
ang12 = [44; 0];
ang13 = [60; 0]

angs = [ang1 ang2 ang3 ang4 ang5 ang6 ang7 ang8 ang9 ang10 ang11 ang12...
    ang13];
ang = [ang1(1,1) ang2(1,1) ang3(1,1) ang4(1,1) ang5(1,1) ang6(1,1)...
    ang7(1,1) ang8(1,1) ang9(1,1) ang10(1,1) ang11(1,1) ang12(1,1)...
    ang13(1,1)]/180*pi;

pos = getElementPosition(ula)/lambda; 

Nsamp = 100;        % 100 samples
nPower = 1/10;       % SNR of 10dB
sources = length(ang);
% Assignment Matrix
A = zeros(sources,11);
for k = 1:sources
    A(k,:) = exp(-1i*2*pi*lambda/2*sin(ang(k))/lambda*[0:10]);
end    
A = transpose(A);
    
rs = rng(1996);     % pseudo random states for reproducible results
signal = sensorsig(pos,Nsamp,angs,nPower); %generated signal model samples

%-----------------------BEAMSCAN-------------------------------------------

spatialspectrum = phased.BeamscanEstimator('SensorArray',ula,...
            'OperatingFrequency',fc,'ScanAngles',-90:90);
spatialspectrum.DOAOutputPort = true;
spatialspectrum.NumSignals = sources;
[~, ang_beamscan] = step(spatialspectrum, signal)
% plotSpectrum(spatialspectrum);
% hold on;
%--------------------MVDR--------------------------------------------------

mvdrspatialspect = phased.MVDREstimator('SensorArray',ula,...
        'OperatingFrequency',fc,'ScanAngles',-90:90,...
        'DOAOutputPort',true,'NumSignals',sources);
[~,ang_mvdr] = step(mvdrspatialspect, signal)
% plotSpectrum(mvdrspatialspect);

%------------------MUSIC---------------------------------------------------

musicspatialspect = phased.MUSICEstimator('SensorArray',ula,...
        'OperatingFrequency',fc,'ScanAngles',-90:90,...
        'DOAOutputPort',true,'NumSignalsSource','Property','NumSignals',min(sources,10));
[~,ang_music] = musicspatialspect(signal)
ymvdr = mvdrspatialspect(signal);
ymusic = musicspatialspect(signal);
helperPlotDOASpectra(mvdrspatialspect.ScanAngles,...
        musicspatialspect.ScanAngles,ymvdr,ymusic,'ULA')
hold all;

%----------------------IG Pencil-------------------------------------------

% signal = transpose(signal);

%Since true covariance matrix is unavailable it is replaced by sample
%covariance matrix
signal_size = size(signal);
R_hat = zeros(signal_size(2),signal_size(2));
for i = 1:signal_size(1)
    R_hat = R_hat + signal(i,:)'*signal(i,:);
end
R_hat = R_hat/Nsamp;
    
% R_hat = signal'*signal;

%eigen values and eigen vectors of sample covariance matrix
[N,V] = eig(R_hat);

%calculating the function values for -90 to 90 with 1 interval for 1 degrees
phi = -90:1:90;
yig = zeros(length(phi),1);
for iter=1:length(phi)
    a_phi = zeros(11,1);
    for j = 1:11
        a_phi(j,1) = exp(-1i*2*pi*(j-1)*uniform_distance/lambda*...
            sin(phi(1,iter)/180*pi));
    end
% As mentioned in the paper, calculation of a_phi'*R-1hat*a_phi
    pencil = abs(a_phi'*inv(R_hat)*a_phi);
    yig(iter,1) = 1/(log(pencil))^2;
end
% calculation of spatial spectrum for IG
spec_ig = yig/max(yig);

spec_ig =10*log10(spec_ig);


% A linear search across the function will give you the position of maximas

[psor,lsor] = findpeaks(yig,phi,'SortStr','descend');
ang_ig = lsor(1:sources)

% Plotting the spectrum of IG
plot(phi, spec_ig, 'k');
plot(ang/pi*180, zeros(length(ang)), 'k*')
legend('MVDR','MUSIC','IG Pencil', 'TRUE DOA');



    
    
        
        
    







