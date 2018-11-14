clear all;
%---------------------SETUP MODEL------------------------------------------
%uniform linear array of 11 sensors with isotropic power kept at 0.5 lambda

c = physconst('LightSpeed');
fc = physconst('LightSpeed');              % Operating frequency
lambda = c/fc;  %wavelength of 1 m
uniform_distance = lambda/2;

ula = phased.ULA('NumElements',11,'ElementSpacing',uniform_distance);

%---------------------signal model-----------------------------------------
% Two sources placed at far locations
ang1 = [-20; 0];          
ang2 = [30; 0];

%Different sets of angles which are closer than the earlier set
ang1 = [-23; 0];
ang2 = [-20; 0];

%signal model with 10 sources at equally spaced between -60 and +60
angs = [ang1 ang2];
ang = [ang1(1,1) ang2(1,1)];

pos = getElementPosition(ula)/lambda; 

Nsamp = 100;        % 100 samples for generated signal 


% We are running monte carlo for 20 discrete values of the SNR 
k = linspace(-20,20,20);
% k = linspace(-5,20,40);
param_power_noise = 10.^(-k/10);

%storing the estimated values for each snr value
estimatedmv1 = zeros(20,1);    
estimatedmu1 = zeros(20,1);    
estimatedigp1 = zeros(20,1);    

estimatedmv2 = zeros(20,1);     
estimatedmu2 = zeros(20,1);     
estimatedigp2 = zeros(20,1);    

% Corresponding to estimated angles, we are also storing the squared sum of
% estimated angles

sq_estimatedmv1 = zeros(20,1);
sq_estimatedmu1 = zeros(20,1);
sq_estimatedigp1 = zeros(20,1);

sq_estimatedmv2 = zeros(20,1);
sq_estimatedmu2 = zeros(20,1);
sq_estimatedigp2 = zeros(20,1);

% Corresponding variances for the estimated angles

variance_mv1 = zeros(20,1);
variance_mu1 = zeros(20,1);
variance_igp1 = zeros(20,1);

variance_mv2 = zeros(20,1);
variance_mu2 = zeros(20,1);
variance_igp2 = zeros(20,1);

% Corresponding minimum mean squared error function for plotting

errormv1 = zeros(20,1);
errormu1 = zeros(20,1);
errorigp1 = zeros(20,1);

errormv2 = zeros(20,1);
errormu2 = zeros(20,1);
errorigp2 = zeros(20,1);

% monte carlo samples
monte_samples = 100;


for i=1:20
    for j=1:monte_samples
        rng(j+100);
        signal = sensorsig(pos,Nsamp,angs,param_power_noise(i)); %generated signal model samples
    
        mvdrspatialspect = phased.MVDREstimator('SensorArray',ula,...
                'OperatingFrequency',fc,'ScanAngles',-90:90,...
                'DOAOutputPort',true,'NumSignals',2);
        [~,ang_mvdr] = step(mvdrspatialspect, signal);
        
        estimatedmv1(i) = estimatedmv1(i) + max(ang_mvdr);
        estimatedmv2(i) = estimatedmv2(i) + min(ang_mvdr);
        
        sq_estimatedmv1(i) = sq_estimatedmv1(i) + max(ang_mvdr)^2;
        sq_estimatedmv2(i) = sq_estimatedmv2(i) + min(ang_mvdr)^2;
        
        errormv1(i) = errormv1(i) + (max(ang_mvdr)-ang(2))^2;
        errormv2(i) = errormv2(i) + (min(ang_mvdr)-ang(1))^2;
        
    %------------------MUSIC---------------------------------------------------

        musicspatialspect = phased.MUSICEstimator('SensorArray',ula,...
                'OperatingFrequency',fc,'ScanAngles',-90:90,...
                'DOAOutputPort',true,'NumSignalsSource','Property','NumSignals',2);
        [~,ang_music] = musicspatialspect(signal);

        estimatedmu1(i) = estimatedmu1(i) + max(ang_music);
        estimatedmu2(i) = estimatedmu2(i) + min(ang_music);
        
        sq_estimatedmu1(i) = sq_estimatedmu1(i) + max(ang_music)^2;
        sq_estimatedmu2(i) = sq_estimatedmu2(i) + min(ang_music)^2;
        
        errormu1(i) = errormu1(i) + (max(ang_music)-ang(2))^2;
        errormu2(i) = errormu2(i) + (min(ang_music)-ang(1))^2;
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

        estimatedigp1(i) = estimatedigp1(i) + max(ang_ig);
        estimatedigp2(i) = estimatedigp2(i) + min(ang_ig);
        
        sq_estimatedigp1(i) = sq_estimatedigp1(i) + max(ang_ig)^2;
        sq_estimatedigp2(i) = sq_estimatedigp2(i) + min(ang_ig)^2;
        
        errorigp1(i) = errorigp1(i) + (max(ang_ig)-ang(2))^2;
        errorigp2(i) = errorigp2(i) + (min(ang_ig)-ang(1))^2;  
        
    end
    estimatedmv1(i) = estimatedmv1(i)/monte_samples;
    estimatedmu1(i) = estimatedmu1(i)/monte_samples;
    estimatedigp1(i) = estimatedigp1(i)/monte_samples;
    estimatedmv2(i) = estimatedmv2(i)/monte_samples;
    estimatedmu2(i) = estimatedmu2(i)/monte_samples;
    estimatedigp2(i) = estimatedigp2(i)/monte_samples;
    
    sq_estimatedmv1(i) = sq_estimatedmv1(i)/monte_samples;
    sq_estimatedmu1(i) = sq_estimatedmu1(i)/monte_samples;
    sq_estimatedigp1(i) = sq_estimatedigp1(i)/monte_samples;
    sq_estimatedmv2(i) = sq_estimatedmv2(i)/monte_samples;
    sq_estimatedmu2(i) = sq_estimatedmu2(i)/monte_samples;
    sq_estimatedigp2(i) = sq_estimatedigp2(i)/monte_samples;
    
    variance_mv1(i) = sq_estimatedmv1(i) - (estimatedmv1(i))^2 ;
    variance_mu1(i) = sq_estimatedmu1(i) - (estimatedmu1(i))^2 ;
    variance_igp1(i) = sq_estimatedigp1(i) - (estimatedigp1(i))^2 ;
    variance_mv2(i) = sq_estimatedmv2(i) - (estimatedmv2(i))^2 ;
    variance_mu2(i) = sq_estimatedmu2(i) - (estimatedmu2(i))^2 ;
    variance_igp2(i) = sq_estimatedigp2(i) - (estimatedigp2(i))^2 ;
    
    errormv1(i) = errormv1(i)/monte_samples;
    errormu1(i) = errormu1(i)/monte_samples;
    errorigp1(i) = errorigp1(i)/monte_samples;
    errormv2(i) = errormv2(i)/monte_samples;
    errormu2(i) = errormu2(i)/monte_samples;
    errorigp2(i) = errorigp2(i)/monte_samples;
   
end    


plot(k,errormv1,'g--+',k,errormu1,'b--o',k,errorigp1,'-.*c',k,errormv2,...
    'r--+',k,errormu2,'k--o',k,errorigp2,'-.*y')
title ('Monte Carlo Comparison of algorithms')
xlabel('SNR (in dB)')
ylabel('MMSE error')
legend('MVDR1','MUSIC1','IGPencil1','MVDR2','MUSIC2','IGPencil2');
