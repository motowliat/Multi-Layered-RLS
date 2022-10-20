clc
clear 
% close all


Cs = 2;% constellation size of mubulation (BPSK)

% double rate
M1 = 50;% must be odd. order of h1, the SI channel
M1_fix = 0;% number of paths of h1 which are fixed and time-invariying 

Max_order = M1;% max of channel orders


lambda = 1-.5/M1; %forgetting factor for fixed lambda
delta = 0.001;% initialer of RLS

Max_Layers = 5;% for multi-layered


B = 5e3;% Bandwidth 
COT1 = 40e-3;% coherence time 
Change_interval = 1;%the channel is changing after this many samples in single rate 

Energy_h1 = 0;% in (dB)


L_train = 500000;% symbol length for trainung
CHANNEL = 20; % number of channels to repeat the process
snr = 4000;% evaluating SNR in (dB)



pdp_h1 = exp(-.4*(0:M1-1));% pdp of SI channel in double rate (w)  


E1 = 10^(Energy_h1/10); %in (w) 
pdp_h1 = E1/sum(pdp_h1) * pdp_h1;% adjusted pdp

noise_power = Energy_h1-snr;% in (dB)
N0 = 10^(noise_power/10); % in (w)

[V1,U1] = AFC_generator(COT1,B,Change_interval);% V1 is the target ACF and U1 the filter to this end

% %-- loading the channel gererator matrices
for channel = 1:CHANNEL
    WW2{channel} = 1/sqrt(2)*(randn(M1,length(U1)+floor((L_train)/Change_interval))+1i*randn(M1,length(U1)+floor((L_train)/Change_interval)));
end
save('WW2','WW2')% Not to interupt with WW1 in the other codes
load('WW2','WW2')



for channel = 1:CHANNEL
    
%%%%%%%%%%%%%% CHANNEL PREPARATION %%%%%%%%%%%%%%%
    W1 = WW2{channel};

    %%%% Training channels generation %%%

    W1_int = 1:length(U1); % interval time of W1 that counts

    hH1_train = zeros(L_train,M1);% h1 channel matrix for training
    for l = 1:L_train-1
        clc
        fprintf('Channel=%g \n Channel preparation... %g of %g \n',channel,l,L_train-1);
        h1_variable = channel_gen(W1(1:M1,W1_int),pdp_h1,U1,M1_fix);
        W1_int = W1_int+1;

        hH1_train(l,:) = h1_variable; % the channel varies just for signle rate not double rate
    end


    %     % Checking the ACF
    V_hat = ACF_estimator(hH1_train(:,1),length(U1));% estimate of ACF
    %         
%         figure
%         plot(V1); hold all
%         plot(abs(V_hat))

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         

     %%               

    X1 = qammod(randi([0,Cs-1],1,L_train),Cs);%

    x1 = [zeros(1,M1-1) , X1];% added zeros symbos to check the BER


    s = zeros(1,L_train-Max_order+1);
    for l = 1:L_train
        h1_variable = hH1_train(l,:);
        s(l) = h1_variable*x1(l+M1-1:-1:l).';%SI signal
    end


    noise = sqrt(N0/2)*(randn(1,L_train-Max_order+1)+1i*randn(1,L_train-Max_order+1)); %ambiant noise

    y = s(1:L_train-Max_order+1)+noise(1:L_train-Max_order+1);% Training with noise


    %--- using the signals to estimate the channels

    H_hat = FixedLambda_MultiLayered3_Checking_ACF(x1,y,M1,delta,lambda,Max_Layers,channel);
   
    
    
    
    
    h{1} = hH1_train(1000:L_train-Max_order+1,:);% original channel
    for c = 1:Max_Layers
        
        h_hat{c} = H_hat(:,:,c);
        
        h{c+1} = h{c}-h_hat{c}(1000:L_train-Max_order+1,:);
        
        for m = 1:10% for the first 10 taps
            V(:,m) = ACF_estimator(h{c+1}(:,m),length(U1)).';% estimate of ACF (for the first tap)
            Test_V1(:,m) = ACF_estimator(h{1}(:,m),length(U1)).'; % For the original channel. just to check if we are doing it right
        end
        VV(:,c+1,channel) = mean(V,2);% average over all cannel taps
        Test_VV1(:,channel) = mean(Test_V1,2);% test
    end
    VV(:,1,channel) = V1(1:length(VV(:,2,1))).';% for the original channel
end

VVV = mean(VV,3);% average over all channels
Test_VVV1 = mean(Test_VV1,2);%test

figure
plot(abs(VVV)); hold all
plot(abs(Test_VVV1),'Linewidth',2);
plot(1:length(VVV(:,1)),0.5*ones(1,length(VVV(:,1))),'r-.')
ylim([0 1.1])
legend('${\varphi _{(1)}}[m]$','${\varphi _{(2)}}[m]$',...
   '${\varphi _{(3)}}[m]$','${\varphi _{(4)}}[m]$',...
   '${\varphi _{(5)}}[m]$','${\varphi _{(6)}}[m]$','${\varphi _{(1)}}[m]$ Test','Threshold','Interpreter','latex')




