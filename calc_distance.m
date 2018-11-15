%Function that takes two matrices as input and returns
% the maximum of the absolute value of the difference between the two.
function d = calc_distance(mu,mu1)

d = max(max(abs(mu-mu1))');

end