function [H1_hat,e,e_f,Lambda,LAYERS] = AdaptiveLambda_MultiLayered3(x1,d,N0,M1,delta,lambda_max,Max_Layers,Input_mode)

Lambda_Adaptivity = 1;
mode = 2;% lambda adaptivity mode.


Layers_Adaptivity = 1;


a = .97; %damping factor

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
           [lambda,sigma_e2,sigma_q2,sigma_v2] = Adaptive_Forget_Factor...
               (sigma_e2,sigma_q2,sigma_v2,N0,e(l),q,M1,lambda_max);
           
        
%         SIGMA_e2(l) = sigma_e2;
%         SIGMA_q2(l) = sigma_q2;
    end
    Lambda(1,l) = lambda;
    
%     % Calculating the a posteriori error
     E(l) = d(l)-w'*uvec; %the error after update
     E_avg = a*E_avg+(1-a)*(abs(E(l))).^2;% obtaing the real-time poer of E(l)
     s(l) = E_avg;% Just to plot

     
%     if l >= 4500
%         l;
%     end
    

    % ----Second layer and after layers, upldate h1
    uvec_f = uvec;
    h1_f = h1;
    h1_hatImp = h1_f;
    r_hat_f(l) = d(l);
    for c = 1:Max_Layers % We need Max_Layers+1 a posteriori error 

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
       if c <= Layers-1% do not count the overfitted results
          h1_hatImp = h1_hatImp + h1_f;% improved estimation of h1   
       end
       
       
%     % Adaptive forgetting factor at RLS#>1
        if Lambda_Adaptivity == 1
           switch mode
               case 1  % (1) Simple mode [C. Paleologu et al.]. lambda is calculated 
                        % for the first layer and it is the same for all
                        % layers
                    lambda_f{c} = lambda;

               case 2  % (2) Complete mode. lambda is different at each layer 
        %                               and it is calculated based on [C.
        %                               Paleologu et al.]              
                    [lambda_f{c},sigma_e2_f{c},sigma_q2_f{c},sigma_v2_f{c}] = Adaptive_Forget_Factor...
                        (sigma_e2_f{c},sigma_q2_f{c},sigma_v2_f{c},N0,e_f{c}(l),q_f{c},M1,lambda_max);
                    Lambda(c+1,l) = lambda_f{c};              
           end
        end
        
        
%     % Calculating the a posteriori error for second layer and after
        E_f{c} = r_hat_f(l)-w_f{c}'*uvec_f; %the error after update at the final layer
        E_avg_f{c} = a*E_avg_f{c}+(1-a)*(abs(E_f{c})).^2;% obtaing the real-time poer of E(l)
        s_f{c}(l) = E_avg_f{c};% Just to plot

    end

    h1 = h1_hatImp;% use the improved h1 estimation via the multilayer, in place of the sigle layer one
    H1_hat(l,:) = h1;   
    
   

% % Adaptive number of layers
    % gathering the a posteriori powers in a sigle vector Pi
     Pi(1) = E_avg;
     for c = 1:Max_Layers
         Pi(c+1) = E_avg_f{c};
     end


% %     Power order checking
%      Layers = 5;
%      if Pi(1) < 1*N0
%         Layers = 1;
%      end

% %   Power order checking
     if Input_mode == 3 % correlated
         Threshold = 1.81*N0;
     else % not correlated
         Threshold = N0;
     end
     Layers = 1;
     for k = 2:Max_Layers+1
         if  (Pi(k-1) > Threshold) && (Pi(k-1) > Pi(k))
            Layers = min(k,Max_Layers);
         else
            break;
         end
     end
     
     
        
     LAYERS(1,l) = Layers;
end




% figure
% hold all
% plot(1*N0*ones(1,length(s)))
% plot(s)
% if Layers_Adaptivity == 1
%     for c = 1:Max_Layers
%         plot(s_f{c})
%     end
% end
% set(gca, 'YScale', 'log')
% xlim([4800,6500])

end














