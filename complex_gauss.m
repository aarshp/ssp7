function p = complex_gauss(z,sigma)

a = [real(z);imag(z)];
sig = [real(sigma),-imag(sigma);imag(sigma),real(sigma)];
p = mvnpdf(a,[0;0;0;0],sig);

end