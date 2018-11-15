function d = calc_distance(phi,phi1,J,K)

d = 0;
% d = max(max(abs(mu-mu1))');

for jj = 1:J
    for k = 1:K
        dd = norm(reshape(phi(jj,k,:,:),[2,2])-reshape(phi1(jj,k,:,:),[2,2]));
        if dd > d
            d = dd;
        end
    end
end

end