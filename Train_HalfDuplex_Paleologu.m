% Layers=1, and Adaptive forgetting factor according to [Paleologu]

function [H1_hat,e,Lambda] = Train_HalfDuplex_Paleologu(x1,d,N0,M1,delta,lambda_max)

Lambda_Adaptivity = 1;


N = length(d);

lambda = lambda_max; % for start


H1_hat = zeros(N,M1);


h1 = zeros(M1,1);
w = conj(h1);% for ZF


% vectors that we use for channel estimation
u1 = x1(:);
d = d(:);


% RLS initialization
p = eye(M1)/delta;
e = zeros(1,N);
y_hat = zeros(1,N);
 


%For adaptive forgetting factor
SIGMA_e2 = zeros(1,N);
SIGMA_q2 = zeros(1,N);
SIGMA_v2 = zeros(1,N);
        
sigma_e2 = 0;
sigma_q2 = 0;
sigma_v2 = 0;


Lambda = zeros(1,N);

for l = 1:N
    
    %-----ZF RLS---------------
        
    uvec = u1(l+M1-1:-1:l);
    k = lambda^(-1)*p*uvec/(1+lambda^(-1)*uvec'*p*uvec);
    y_hat(l) = w'*uvec;
    e(l) = d(l)-y_hat(l);
    w = w+k*conj(e(l));
    q = uvec'*p*uvec;% for adaptive lambda
    p = lambda^(-1)*p-lambda^(-1)*k*uvec'*p;
        
    h = conj(w);% channel is the conjugate of w
    h1 = h;
    
%    % Adaptive forgetting factor at RLS#1
    if Lambda_Adaptivity == 1
       [lambda,sigma_e2,sigma_q2,sigma_v2] = Adaptive_Forget_Factor_Paleologu...
               (sigma_e2,sigma_q2,sigma_v2,e(l),q,N0,M1,lambda_max);
           
        Lambda(1,l) = lambda;
        SIGMA_e2(l) = sigma_e2;
        SIGMA_q2(l) = sigma_q2;
        SIGMA_v2(l) = sigma_v2;
    end
   
   H1_hat(l,:) = h1;
   
   
       if l > 3000
            l;
       end 
end




end














