clear all;
%---------------------SETUP MODEL------------------------------------------
%uniform linear array of 11 sensors with isotropic power kept at 0.5 lambda

c = physconst('LightSpeed');
fc = physconst('LightSpeed');              % Operating frequency
lambda = c/fc;  %wavelength of 1 m
uniform_distance = lambda/2;

ula = phased.ULA('NumElements',11,'ElementSpacing',uniform_distance);

%---------------------signal model-----------------------------------------
ang1 = [-20; 0];          % First signal
ang2 = [30; 0];

%signal model with 10 sources at equally spaced between -60 and +60
angs = [ang1 ang2];
ang = [ang1(1,1) ang2(1,1)]/180*pi;

pos = getElementPosition(ula)/lambda; 

Nsamp = 100;        % 100 samples
k = linspace(-20,20,20);
param_power_noise = 10.^(-linspace(-20, 20, 20)/10);

errormv1 = zeros(20,1);
errormu1 = zeros(20,1);
errorigp1 = zeros(20,1);

errormv2 = zeros(20,1);
errormu2 = zeros(20,1);
errorigp2 = zeros(20,1);

for i=1:20
    for j=1:10
        rng(j);
        signal = sensorsig(pos,Nsamp,angs,param_power_noise(i)); %generated signal model samples
    
        mvdrspatialspect = phased.MVDREstimator('SensorArray',ula,...
                'OperatingFrequency',fc,'ScanAngles',-90:90,...
                'DOAOutputPort',true,'NumSignals',2);
        [~,ang_mvdr] = step(mvdrspatialspect, signal);
    % plotSpectrum(mvdrspatialspect);
        errormv1(i) = errormv1(i) + (max(ang_mvdr)-30)^2;
        errormv2(i) = errormv2(i) + (min(ang_mvdr)+20)^2;
        
    %------------------MUSIC---------------------------------------------------

        musicspatialspect = phased.MUSICEstimator('SensorArray',ula,...
                'OperatingFrequency',fc,'ScanAngles',-90:90,...
                'DOAOutputPort',true,'NumSignalsSource','Property','NumSignals',2);
        [~,ang_music] = musicspatialspect(signal);
        errormu1(i) = errormu1(i) + (max(ang_music)-30)^2;
        errormu2(i) = errormu2(i) + (min(ang_music)+20)^2;
    %----------------------IG Pencil-------------------------------------------

    % signal = transpose(signal);

    %Since true covariance matrix is unavailable it is replaced by sample
    %covariance matrix
        R_hat = signal'*signal;
    
        %eigen values and eigen vectors of sample covariance matrix
        [N,V] = eig(R_hat);
    
        %Linear Search for maximas from -pi/2 to pi/2
        phi = -90:1:90;
        yig = zeros(length(phi),1);
        for iter=1:length(phi)
            a_phi = zeros(11,1);
            for j = 1:11
                a_phi(j,1) = exp(-1i*2*pi*(j-1)*uniform_distance/lambda*...
                    sin(phi(1,iter)/180*pi));
            end
            pencil = abs(a_phi'*inv(R_hat)*a_phi);
            yig(iter,1) = 1/pencil;
        end
        [psor,lsor] = findpeaks(yig,phi,'SortStr','descend');
        ang_ig = lsor(1:2);
        errorigp1(i) = errorigp1(i) + (max(ang_ig)-30)^2;
        errorigp2(i) = errorigp2(i) + (min(ang_ig)+20)^2;
    end
    errormv1(i) = errormv1(i)/1000;
    errormu1(i) = errormu1(i)/1000;
    errorigp1(i) = errorigp1(i)/1000;
    errormv2(i) = errormv2(i)/1000;
    errormu2(i) = errormu2(i)/1000;
    errorigp2(i) = errorigp2(i)/1000;
end    
plot(k,errormv1,'g',k,errormu1,'b--o',k,errorigp1,'-.*c',k,errormv2,'r',k,errormu2,'k--o',k,errorigp2,'-.*y')
legend('MVDR1','MUSIC1','IGPencil1','MVDR2','MUSIC2','IGPencil2');
