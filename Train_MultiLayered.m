% Layers=1, and Adaptive forgetting factor

function [H1_hat,e,e_f,Lambda,LAYERS] = Train_MultiLayered(x1,d,N0,M1,delta,lambda_max,Layers_Adaptivity,Max_Layers)


Lambda_Adaptivity = 0;
mode = 1;% lambda adaptivity mode.


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

% Adaptive number of layers
E = zeros(1,N);
E_avg = 0;


r_hat_f = zeros(N,1);   

% for the second layer and after (in the mulltilayer estimator, the first
% layer is performed separately in this function
for c = 1:Max_Layers
    w_f{c} = zeros(M1,1);
    e_f{c} = zeros(M1,1);
    k_f{c} = zeros(M1,1);
    p_f{c} = eye(M1)/delta;
    
    E_avg_f{c} = 0;
    sigma_e2_f{c} = 0;
    sigma_q2_f{c} = 0;
    sigma_v2_f{c} = 0;
    sigma_b2_f{c} = 0;
    lambda_f{c} = lambda;
end

if Layers_Adaptivity == 1
    Layers = 1;% at the begining the number of layers is 1
else
    Layers = Max_Layers;% if the adaptivity is off, we take the Max_Layers as the number of layers
end

%For adaptive forgetting factor
sigma_e2 = 0;
sigma_q2 = 0;
sigma_v2 = 0;



Lambda = zeros(Max_Layers,N);
LAYERS = zeros(1,N);
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
           [lambda,sigma_e2,sigma_q2,sigma_v2,sigma_b2] = Adaptive_Forget_Factor...
               (sigma_e2,sigma_q2,sigma_v2,sigma_b2,N0,e(l),q,M1,lambda_max);
           
        Lambda(1,l) = lambda;
        SIGMA_e2(l) = sigma_e2;
        SIGMA_q2(l) = sigma_q2;
        SIGMA_b2(l) = sigma_b2;
    end

    
%     % Adaptive number of layers
     E(l) = d(l)-w'*uvec; %the error after update
     E_avg = .995*E_avg+(1-.995)*(abs(E(l))).^2;% obtaing the real-time poer of E(l)
     s(l) = E_avg;% Just to plot ?????????????????
     if (E_avg > N0) && (Layers_Adaptivity == 1)
         Layers = max(2,Layers);
     elseif  Layers_Adaptivity == 1
         Layers = 1;
     end
     
     
% if l >= 5000
%     l;
% end
    

    % ----Second layer and after layers, upldate h1
    uvec_f = uvec;
    h1_f = h1;
    h1_hatImp = h1_f;
    r_hat_f(l) = d(l);
    for c = 1:Layers-1 % The second layer and after. The first SI cancellation is performed above, separately

       s_hat_f = h1_f.'*u1(l+M1-1:-1:l); % the residual of s_hat(l) at the c_th iteration
       r_hat_f(l) = r_hat_f(l)-s_hat_f;

       % RLS
       k_f{c} = lambda_f{c}^(-1)*p_f{c}*uvec_f/(1+lambda_f{c}^(-1)*uvec_f'*p_f{c}*uvec_f);
       y_tild_f = w_f{c}'*uvec_f;
       e_f{c}(l) = r_hat_f(l)-y_tild_f;
       w_f{c} = w_f{c}+k_f{c}*conj(e_f{c}(l)); % ^^^^ WE NEED THIS FOR ADATIVE LAMBDA ^^^^
       q_f{c} = uvec'*p_f{c}*uvec;% for adaptive lambda
       p_f{c} = lambda_f{c}^(-1)*p_f{c}-lambda_f{c}^(-1)*k_f{c}*uvec_f'*p_f{c};


       h_f = conj(w_f{c});
       h1_f = h_f(1:M1);     
       h1_hatImp = h1_hatImp + h1_f;% improved estimation of h1   
       
       
       
%     % Adaptive forgetting factor at RLS#>1
        if Lambda_Adaptivity == 1
           switch mode
               case 1  % (1) Simple mode [C. Paleologu et al.]. lambda for all layers
    %                is the same and it is calculated based on the paper b{1}=0.
                    lambda_f{c} = lambda;

               case 3  % (3) Complete mode (our derivation). lambda is different at each layer 
        %                               and it is calculated based on [C.
        %                               Paleologu et al.] (b=0 in all
        %                               cases)               
                    [lambda_f{c},sigma_e2_f{c},sigma_q2_f{c},sigma_v2_f{c},sigma_b2_f{c}] = Adaptive_Forget_Factor...
                        (sigma_e2_f{c},sigma_q2_f{c},sigma_v2_f{c},sigma_b2_f{c},N0,e_f{c}(l),q_f{c},M1,lambda_max);
                    Lambda(c+1,l) = lambda_f{c};              
           end
        end
    
    
       
%     % Adaptive number of layers
        E_f{c} = r_hat_f(l)-w_f{c}'*uvec_f; %the error after update
        E_avg_f{c} = .995*E_avg_f{c}+(1-.995)*(abs(E_f{c})).^2;% obtaing the real-time poer of E(l)
        s_f{c}(l) = E_avg_f{c};% Just to plot%?????????????
        if (E_avg_f{c} > N0) && (Layers_Adaptivity == 1)
           Layers = min(c+2,Max_Layers);
        end
    
     
    end

    h1 = h1_hatImp;% use the improved h1 estimation via the multilayer, in place of the sigle layer one

    H1_hat(l,:) = h1;   
    
   
LAYERS(1,l) = Layers;
end



end














