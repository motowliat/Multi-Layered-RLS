function [H1_hat,e,Lambda] = PaleologuAdaptiveLambda_RLS(x1,d,N0,M1,delta,lambda_max)


Lambda_Adaptivity = 1;

a = .96; %damping factor

N = length(d);

lambda = lambda_max;% for start

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

% a posteriori error initiation
E = zeros(1,N);
E_avg = 0;




%For adaptive forgetting factor
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
    p = lambda^(-1)*p-lambda^(-1)*k*uvec'*p;
        
    h = conj(w);% channel is the conjugate of w
    h1 = h;
    H1_hat(l,:) = h1; 

%    % Adaptive forgetting factor at RLS#1
    if Lambda_Adaptivity == 1
        if l > 500 % use the Paleologu method after the initial error is setteled
           q = uvec'*p*uvec;% for adaptive lambda
               [lambda,sigma_e2,sigma_q2,sigma_v2] = Adaptive_Forget_Factor_Paleologu...
                   (sigma_e2,sigma_q2,sigma_v2,e(l),q,N0,M1,lambda_max);


            SIGMA_e2(l) = sigma_e2;
            SIGMA_q2(l) = sigma_q2;
            SIGMA_v2(l) = sigma_v2;
        end
    end
    Lambda(1,l) = lambda;
    
    %  Calculating the a posteriori error
     E(l) = d(l)-w'*uvec; %the error after update
     E_avg = a*E_avg+(1-a)*(abs(E(l))).^2;% obtaing the real-time poer of E(l)
     s(l) = E_avg;% Just to plot 



%     if l >= 1101
%         l;
%     end

end




% figure
% % subplot(3,1,1)
% hold all
% plot(N0*ones(1,length(s)))
% plot(s)
% set(gca, 'YScale', 'log')
% xlim([4500,6500])

end

