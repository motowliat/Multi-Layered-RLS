% Layers=1, and Adaptive forgetting factor

function [H1_hat,e] = NLMS(x1,d,M1,miu)

% miu: The LMS step size 

a = .96; %damping factor

N = length(d);

H1_hat = zeros(N,M1);

h1 = zeros(M1,1);
w = conj(h1);% for ZF


% vectors that we use for channel estimation
u1 = x1(:);
d = d(:);


% LMS initialization
e = zeros(1,N);
y_hat = zeros(1,N);

% a posteriori error initiation
E = zeros(1,N);
E_avg = 0;  

for l = 1:N
    
    %-----LMS--------------
    uvec = u1(l+M1-1:-1:l);
    y_hat(l) = w'*uvec;
    e(l) = d(l)-y_hat(l);
    w = w+miu/(uvec'*uvec)*uvec*conj(e(l));% normalized LMS
        
    h = conj(w);% channel is the conjugate of w
    h1 = h;
    H1_hat(l,:) = h1; 


    %  Calculating the a posteriori error
     E(l) = d(l)-w'*uvec; %the error after update
     E_avg = a*E_avg+(1-a)*(abs(E(l))).^2;% obtaing the real-time poer of E(l)
     s(l) = E_avg;% Just to plot 



    if l >= 4500
        l;
    end

end




% figure
% % subplot(3,1,1)
% hold all
% plot(N0*ones(1,length(s)))
% plot(s)
% set(gca, 'YScale', 'log')
% xlim([4500,6500])

end














