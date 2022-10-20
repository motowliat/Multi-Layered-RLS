function [lambda,sigma_k] = Adaptive_Forget_Factor_Bhotto(e_k,sigma_k,N0,M)

b = .99;
c1 = 4;
sigma_v = sqrt(N0);
t = sqrt(c1*sigma_v^2);

e_kf = sign(real(e_k))*max(abs(e_k)-t,0);
sigma_k = b*sigma_k+(1-b)*abs(e_kf)^2;
lambda = 1-2*sigma_k/(M*(sigma_k+c1*sigma_v^2));