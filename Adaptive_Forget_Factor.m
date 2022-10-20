%adaptive forgetting factor for RLS algorithm

function [lambda,sigma_e2,sigma_q2,sigma_v2] = Adaptive_Forget_Factor(sigma_e2,sigma_q2,sigma_v2,N0,e,q,L,lambda_max)

K_a = 2;
K_B = 5*K_a;
gamma = 1.5;
epsilon = 1e-8;

alpha = 1-1/(K_a*L);%.98;%
betta = 1-1/(K_B*L);%.995;%

sigma_e2 = alpha*sigma_e2+(1-alpha)*(abs(e)^2);
sigma_q2 = alpha*sigma_q2+(1-alpha)*(abs(q)^2);
sigma_v2 = betta*sigma_v2+(1-betta)*(abs(e)^2);%


sigma_e = sqrt(sigma_e2);
sigma_q = sqrt(sigma_q2);
sigma_v = sqrt(sigma_v2);

if sigma_e <= gamma*sigma_v
    lambda = lambda_max;
else
    lambda = min((sigma_q*sigma_v) / (epsilon+abs(sigma_e-sigma_v)) , lambda_max);
    lambda = max(lambda,.9);% lower boundary for lambda
end


