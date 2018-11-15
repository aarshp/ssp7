function prob = SSP_EM(mic_pos,particle_pos,z,epsilon)

% t : a particular time instant we are working on
% mic_pos =2X3
% particle_pos = 4 x J

c = 340;
% z = 2xK

Ts = 1/8000;
K = size(z);
K = K(2);                    %K : Frequency bins count

J = size(particle_pos,2);
prob = ones(1,J);         % J : Sampled positions for source

    
% h and tau are given
tau = zeros(K,2,2);
for k = 1:K
    xx = sinc(2*pi*k*norm(mic_pos(1,:)-mic_pos(2,:))/(K*Ts*c));
    for ii = 1:2
         for jj = 1:2
             tau(k,ii,jj) = xx;
             if ii==jj
                tau(k,ii,jj) = xx + epsilon;
             end
         end
    end
end

h = zeros(J,K,2);
for jj = 1:J
    for k = 1:K
        for ii = 1:2
            cc = mic_pos(1,:)/2+mic_pos(2,:)/2;
            gamma = atan(particle_pos(4,jj)/particle_pos(3,jj))-atan((mic_pos(ii,2)-particle_pos(2,jj))/(mic_pos(ii,1)-particle_pos(1,jj)));
            h(jj,k,ii) = exp(2j*pi*k*norm(mic_pos(ii,:)-cc)*cos(gamma)/(K*Ts*c));   % assuming centre is the reference microphone then gamma_m = 0 or 180
        end
    end
end


shift = 10000;  % a big number
iter = 0;       % counter
epsil = 0.000001; % percision
max_iter = 50;

phi_y1 = zeros(J,K);
phi1 = zeros(J,K,2,2);     % _init_ phi, sai, phi_r and phi_y as zeros tensors
phi_r1 = zeros(J,K);   % phi_y and phi_r are 1X1 matrices; phi(t,k) is a 2X2
sai = ones(J,1)/J;        % for each j and k
b = rand(2,J);

for k = 1:K
    for jj = 1:J
        b(:,jj) = inv(reshape(tau(k,:,:),[2,2]))*reshape(h(jj,k,:),[2,1])/(reshape(h(jj,k,:),[2,1])'*inv(reshape(tau(k,:,:),[2,2]))*reshape(h(jj,k,:),[2,1]));
        phi_r1(jj, k) = real((z(:,k))'*(eye(2) - b(:,jj)*reshape(h(jj,k,:),[2,1])')*inv(reshape(tau(k,:,:),[2,2]))*(z(:,k)));
        phi_y1(jj, k) = real(((b(:,jj))'*((z(:,k))*z(:,k)'-real(phi_r1(jj,k))*reshape(tau(k,:,:),[2,2]))*b(:,jj)));
        phi1(jj,k,:,:) = reshape(h(jj,k,:),[2,1])*reshape(h(jj,k,:),[2,1])'*phi_y1(jj,k) + reshape(tau(k,:,:),[2,2])*phi_r1(jj,k,:,:);
        r = rank(reshape(phi1(jj,k,:,:),[2,2]));
        if r/2<1
            disp('aa rha hai yaha')
            disp(r)
        end
    end
end
                   
%formatSpec = 'iteration: %d, error: %2.4f, mu1: [%2.4f %2.4f], mu2: [%2.4f %2.4f] \n';

mu = ones(J,K)/J;              % mu : is an J X K dimen matrix    
 
while shift > epsil && iter  < max_iter
    %reshape(phi(1,1,:,:),[2,2])
    %tic;
    mu_new = ones(J,K)/J;
    iter = iter + 1;
    %disp(iter)
    if iter > 0
        for k = 1:K
            for jj = 1:J  
                try 
                    p = complex_gauss(z(:,k),reshape(phi1(jj,k,:,:),[2,2]));
                catch err
                    disp(err);
                    disp(reshape(phi1(jj,k,:,:),[2,2]));
                    disp(jj);
                    disp(k);
                    disp("J,k values");
                end    
                mu_new(jj,k) = sai(jj,1)*p;    % completed it with reshape(phi(jj,k,:,:),[2,2])
            end
            mu_new(:,k) = mu_new(:,k)/sum(mu_new(:,k));
        end
    end
    

    
    for jj = 1:J
        sai(jj,1) = sum(mu_new(jj,:))/K;
    end
    
    if iter>1
        shift = calc_distance(mu,mu_new);
        disp(shift)
        mu = mu_new;
    end
    toc
end

for jj = 1:J                %confirm
    prob(1,jj) = prob(1,jj)*sai(jj,1);
    for k = 1:K
        prob(1,jj) = prob(1,jj)*complex_gauss(z(:,k),reshape(phi1(jj,k,:,:),[2,2]));
    end
end

end

