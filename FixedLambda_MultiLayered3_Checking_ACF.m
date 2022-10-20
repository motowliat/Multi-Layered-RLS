function [HH] = FixedLambda_MultiLayered3_Checking_ACF(x1,d,M1,delta,lambda,Max_Layers,channel)

Layers = Max_Layers;


N = length(d);



h1 = zeros(M1,1);
w = conj(h1);% for ZF


% vectors that we use for channel estimation
u1 = x1(:);
d = d(:);


% RLS initialization
p = eye(M1)/delta;

HH = zeros(N,M1,Max_Layers);

 

% for the second layer and after (in the mulltilayer estimator, the first
% layer is performed separately in this function
for c = 1:Max_Layers
    w_f = zeros(M1,Max_Layers);
    k_f = zeros(M1,Max_Layers);
    p_f(:,:,c) = eye(M1)/delta;
    
end




for l = 1:N
    clc
    fprintf('ACF estimation... %g of %g | channel=%g \n',l,N,channel);
    %-----ZF RLS---------------
        
    uvec = u1(l+M1-1:-1:l);
    k = lambda^(-1)*p*uvec/(1+lambda^(-1)*uvec'*p*uvec);
    y_hat = w'*uvec;
    e = d(l)-y_hat;
    w = w+k*conj(e);
    p = lambda^(-1)*p-lambda^(-1)*k*uvec'*p;
        
    HH(l,:,1) = w';%SAVE THE CHANNELS IN A MATRIX


   

    % ----Second layer and after layers, upldate h1
    uvec_f = uvec;
    h1_f = h1;
    h1_hatImp = h1_f;
    r_hat_f = d(l);
    for c = 1:Max_Layers % We need Max_Layers+1 a posteriori error 

       s_hat_f = h1_f.'*u1(l+M1-1:-1:l); % the residual of s_hat(l) at the c_th iteration
       r_hat_f = r_hat_f-s_hat_f;

       % RLS
       k_f(:,c) = lambda^(-1)*p_f(:,:,c)*uvec_f/(1+lambda^(-1)*uvec_f'*p_f(:,:,c)*uvec_f);
       y_tild_f = w_f(:,c)'*uvec_f;
       e_f = r_hat_f-y_tild_f;
       w_f(:,c) = w_f(:,c)+k_f(:,c)*conj(e_f); % ^^^^ WE NEED THIS FOR ADATIVE LAMBDA ^^^^
       p_f(:,:,c) = lambda^(-1)*p_f(:,:,c)-lambda^(-1)*k_f(:,c)*uvec_f'*p_f(:,:,c);


       h_f = conj(w_f(:,c));
       h1_f = h_f(1:M1); 
       
       if c <= Layers-1% do not count the overfitted results
          h1_hatImp = h1_hatImp + h1_f;% improved estimation of h1 
  
          HH(l,:,c) = h1_f.'; %SAVE THE CHANNELS IN A MATRIX
       end
       
       
    end
    
   

end

end

