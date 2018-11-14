clear;clc;
%Source Motion Model Parameters
delta_t = 0.375;
N_particles = 1000;
N_steps = int64(10/delta_t);
beta = 2;
v_bar = 1;
a = exp(-beta*delta_t);
b = v_bar*sqrt(1-a^2);
F = [eye(2), a*delta_t*eye(2);zeros(2,2),a*eye(2)];
Q = [b^2*delta_t^2*eye(2),zeros(2,2);zeros(2,2),b^2*eye(2)];



%True source and microphone motion models
source = zeros(N_steps+1,4);
source(1,1:2) = 0.5*[randi([3 9],1,2)] ;
mic = zeros(N_steps+1,2,3);
mic(1,1,:) = [1.35,1,1.5];
mic(1,2,:) = [1.65,1,1.5];
mic(end,1,:)  = [4.85,1,1.5];
mic(end,2,:) = [5.15,1,1.5];

for i = 1:3
	mic(:,1,i) = linspace(mic(1,1,i),mic(end,1,i),N_steps+1);
	mic(:,2,i) = linspace(mic(1,2,i),mic(end,2,i),N_steps+1);
end

%rir parameters and model
room_dim = [6 6 2.5];
fs = 8000;
c = 340;
rt60  = 0.5;
rir_samples = fs*delta_t;
%rir(1,:,:) = rir_generator(c,fs,mic(1,:,:),[source(1,1:2) 0],room_dim,rt60,rir_samples);

%Sampled Source Positions
source_samp = zeros(N_steps+1,N_particles,4);
source_samp(1,:,:) = [0.5*randi([3 9],N_particles,2) zeros(N_particles,2)];


%Grid points
[X,Y] = meshgrid(0:0.1:6,0:0.1:6);
grid_pts = [X(:),Y(:)];
kdeprob = zeros(N_steps,size(grid_pts,1));

%Signal Array : sig 
signal_raw = load('timit_audio.mat');
signal_raw = signal_raw.audio_samps;
%signal = reshape(signal_raw(1:(N_steps-1)*fs*delta_t), [N_steps-1 fs*delta_t]);
%signal(N_steps,:) = reshape(signal_raw((N_steps-1)*fs*delta_t +1 :end),[1 (fs*10 -(N_steps-1)*fs*delta_t)]);

%STFT params
window_size = 256;
hop = window_size/4;
window = rectwin(window_size);
n_bins = 2^nextpow2(fs*delta_t);
outside_source  = 0;
outside_samples = 0;


for t = 1:N_steps
	temp_source = mvnrnd(F*reshape(source(t,:),[4 1]),Q,1);
    a = temp_source(1:2) <0 ;
    b = temp_source(1:2) > 6;
    outside_source = outside_source + max(sum(a(:)),sum(b(:)));
	%while ((sum(a(:)) > 0) || (sum(b(:)) > 0))
	%	temp_source = mvnrnd(F*reshape(source(t,:),[4 1])   ,Q,1);
    %end
    
	source(t+1,:) = temp_source;
    if t == N_steps
        rir = rir_generator(c,fs,reshape(mic(t,:,:),[2,3]),[source(t,1:2) 0],room_dim,rt60,fs*10 - fs*(N_steps-1)*delta_t);
        S1 = conv(rir(1,:),signal_raw(fs*(N_steps-1)*delta_t+1:end));
        S2 = conv(rir(2,:),signal_raw(fs*(N_steps-1)*delta_t+1:end));
    else    
        rir = rir_generator(c,fs,reshape(mic(t,:,:),[2,3]),[source(t,1:2) 0],room_dim,rt60,rir_samples);	
        S1 = conv(rir(1,:),signal_raw((t-1)*fs*delta_t + 1: t*fs*delta_t));
        S2 = conv(rir(2,:),signal_raw((t-1)*fs*delta_t + 1: t*fs*delta_t));
    end
    
    
    %STFT with signal array
	%Z1 and 
    [Z1,farr1,tarr1] = stft(S1,window,hop,n_bins,fs);
    [Z2,farr2,tarr2] = stft(S2,window,hop,n_bins,fs);
    
	for j = 1:N_particles
		temp_samp = mvnrnd(F*reshape(source_samp(t,j,:),[4,1]),Q,1);
        a = temp_samp(:,1:2) < 0;
        b =  temp_samp(:,1:2) > 6;
        outside_samples = outside_samples + max(sum(a(:)),sum(b(:)));
		%while ( (sum(a(:)) > 0) || (sum(b(:)) > 0))
		%	temp_samp = mvnrnd(F*reshape(source_samp(t,j,:),[4,1]),Q,1);
        %end
   		source_samp(t+1,j,:) = temp_samp;
	%Update equations for weight
	%w(t+1,j,:) = w(t,j,:)
    end
    
	%kdeprob(i,:),pts = ksdensity(reshape(source_samp(t,:,1:2),[N_samples,2]),grid_pts,'Weights',reshape(w(t,:,:),[N_samples,2]));
end
disp(outside_source);
disp(outside_samples);

%contour(X,Y,reshape(kdeprob(1,:),size(X)));
		
