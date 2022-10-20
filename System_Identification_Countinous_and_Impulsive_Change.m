clc
clear 
% close all


Cs = 2;% constellation size of mubulation (BPSK)

% double rate
M1 = 50;% must be odd. order of h1, the SI channel
M1_fix = 0;% number of paths of h1 which are fixed and time-invariying 

Max_order = M1;% max of channel orders


lambda_max = 1-.1/M1; % max forgetting factor for adaptive lambda
lambda = 1-.5/M1; %forgetting factor for fixed lambda
delta = 0.001;% initialer of RLS

Max_Layers = 5;% for multi-layered


Input_mode = 1;%(1) BPSK, (2) AWGN, (3) AR Process


B = 5e3;% Bandwidth 
COT1 = 40e-3;% coherence time for non-direct paths of the SI channel
Change_interval = 1;%the channel is changing after this many samples in single rate 

Energy_h1 = 0;% in (dB)


L_train = 3000;% symbol length for trainung
CHANNEL = 5; % number of channels to repeat the process
SNR = 10;%[10 20 30];% evaluating RT-SNR in (dB)
ITR = 1; % number of repeating the process for a fixed snr and channel







% pdp_h1 = exp(-.4*(0:M1-1));% pdp of SI channel in double rate (w)  
load('pdp_h1_LakeExperiment'); % pdp of SI channel in double rate (w)  
pdp_h1 = pdp_h1(1:M1); 


E1 = 10^(Energy_h1/10); %in (w) 
pdp_h1 = E1/sum(pdp_h1) * pdp_h1;% adjusted pdp


[V1,U1] = AFC_generator(COT1,B,Change_interval);% autocorrelation function of non-direct paths, SI channel

% %-- loading the channel gererator matrices
% for channel = 1:CHANNEL
%     WW1{channel} = 1/sqrt(2)*(randn(M1,length(U1)+floor((L_train)/Change_interval))+1i*randn(M1,length(U1)+floor((L_train)/Change_interval)));
% end
% save('WW1','WW1')
load('WW1','WW1')


Channel_Est_Error2_0 = zeros(CHANNEL,length(SNR),L_train-Max_order+1);
Channel_Est_Error2_1 = zeros(CHANNEL,length(SNR),L_train-Max_order+1);
Channel_Est_Error2_2 = zeros(CHANNEL,length(SNR),L_train-Max_order+1);
Channel_Est_Error2_3 = zeros(CHANNEL,length(SNR),L_train-Max_order+1);
Channel_Est_Error2_33 = zeros(CHANNEL,length(SNR),L_train-Max_order+1);
Channel_Est_Error2_4 = zeros(CHANNEL,length(SNR),L_train-Max_order+1);
Channel_Est_Error2_5 = zeros(CHANNEL,length(SNR),L_train-Max_order+1);
for channel = 1:CHANNEL
    
    
    
    %%%%%%%%%%%%%% CHANNEL PREPARATION %%%%%%%%%%%%%%%
    W1 = WW1{channel};
     
    %%%% Training channels generation %%%
      
    W1_int = 1:length(U1); % interval time of W1 that counts
    
    hH1_train = zeros(L_train,M1);% h1 channel matrix for training
    for l = 1:L_train-1
        if l == 1 || rem(l,Change_interval) == 0 % because l is the double rate counter
           h1_variable = channel_gen(W1(1:M1,W1_int),pdp_h1,U1,M1_fix);
           W1_int = W1_int+1;
        end
        if l < 1000 % the impulsive change happens at l=500
           hH1_train(l,:) = h1_variable; % the channel varies just for signle rate not double rate
        else
           hH1_train(l,:) = -h1_variable;
        end
        
    end   
                 

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
     
 
 

 
 %%           
         
    
    snr_ctr = 0;
    for snr = SNR
        snr_ctr = snr_ctr+1;
        noise_power = Energy_h1-snr;% in (dB)
        N0 = 10^(noise_power/10); % in (w)
        
        Channel_Est_Error1_0 = zeros(ITR,L_train-Max_order+1);
        Channel_Est_Error1_1 = zeros(ITR,L_train-Max_order+1);
        Channel_Est_Error1_2 = zeros(ITR,L_train-Max_order+1);
        Channel_Est_Error1_3 = zeros(ITR,L_train-Max_order+1);
        Channel_Est_Error1_33 = zeros(ITR,L_train-Max_order+1);
        Channel_Est_Error1_4 = zeros(ITR,L_train-Max_order+1);
        Channel_Est_Error1_5 = zeros(ITR,L_train-Max_order+1);
        for itr = 1:ITR

            
        %%%% Training Phase
            switch Input_mode 
                case 1 % BPSK
                    X1 = qammod(randi([0,Cs-1],1,L_train),Cs);%
                case 2 % AWGN
                    X1 = 1/sqrt(2)*(randn(1,L_train)+1i*randn(1,L_train));
                case 3 % AR Process
                    awgn = 1/sqrt(2)*(randn(1,L_train)+1i*randn(1,L_train));
                    X1 = filter(1,[1,-.9],awgn);
                    X1 = X1./sqrt(mean((abs(X1)).^2));
             end

            x1 = [zeros(1,M1-1) , X1];% added zeros symbos to check the BER

            
            s = zeros(1,L_train-Max_order+1);
            for l = 1:L_train
                h1_variable = hH1_train(l,:);
                s(l) = h1_variable*x1(l+M1-1:-1:l).';%SI signal
            end

            
            noise = sqrt(N0/2)*(randn(1,L_train-Max_order+1)+1i*randn(1,L_train-Max_order+1)); %ambiant noise
            
            y = s(1:L_train-Max_order+1)+noise(1:L_train-Max_order+1);% Training with noise

            
            %--- using the signals to train the channels
            [H1_hat_0,e_0,e_f_0,Lambda_0,LAYERS_0] = BhottoAdaptiveLambda_MultiLayered(x1,y,N0,M1,delta,lambda_max,Max_Layers,Input_mode);
            [H1_hat_1,e_1,e_f_1,Lambda_1,LAYERS_1] = PaleologuAdaptiveLambda_MultiLayered(x1,y,N0,M1,delta,lambda_max,Max_Layers,Input_mode);
            [H1_hat_2,e_2,e_f_2,Lambda_2,LAYERS_2] = FixedLambda_MultiLayered3(x1,y,N0,M1,delta,lambda,Max_Layers,Input_mode);
            [H1_hat_3,e_3,Lambda_3] = BhottoAdaptiveLambda_RLS(x1,y,N0,M1,delta,lambda_max); % simple RLS
            [H1_hat_33,e_33,Lambda_33] = PaleologuAdaptiveLambda_RLS(x1,y,N0,M1,delta,lambda_max); % simple RLS
            [H1_hat_4,e_4,Lambda_4] = FixedLambda_RLS(x1,y,N0,M1,delta,lambda); % simple RLS
            [H1_hat_5,e_5] = NLMS(x1,y,M1,.009*M1);% LMS

            Channel_Est_Error1_0(itr,:) = sum(((abs(hH1_train(1:L_train-Max_order+1,:)-H1_hat_0(1:L_train-Max_order+1,:))).^2).')./E1;
            Channel_Est_Error1_1(itr,:) = sum(((abs(hH1_train(1:L_train-Max_order+1,:)-H1_hat_1(1:L_train-Max_order+1,:))).^2).')./E1;
            Channel_Est_Error1_2(itr,:) = sum(((abs(hH1_train(1:L_train-Max_order+1,:)-H1_hat_2(1:L_train-Max_order+1,:))).^2).')./E1;
            Channel_Est_Error1_3(itr,:) = sum(((abs(hH1_train(1:L_train-Max_order+1,:)-H1_hat_3(1:L_train-Max_order+1,:))).^2).')./E1;
            Channel_Est_Error1_33(itr,:) = sum(((abs(hH1_train(1:L_train-Max_order+1,:)-H1_hat_33(1:L_train-Max_order+1,:))).^2).')./E1;
            Channel_Est_Error1_4(itr,:) = sum(((abs(hH1_train(1:L_train-Max_order+1,:)-H1_hat_4(1:L_train-Max_order+1,:))).^2).')./E1;
            Channel_Est_Error1_5(itr,:) = sum(((abs(hH1_train(1:L_train-Max_order+1,:)-H1_hat_5(1:L_train-Max_order+1,:))).^2).')./E1;
 
            LA0_1(itr,:) = LAYERS_0;
            LA1_1(itr,:) = LAYERS_1;
            LA2_1(itr,:) = LAYERS_2;

            clc
            fprintf('itr=%g , snr=%g , channel=%g\n',itr,snr,channel);

                
        end

        Channel_Est_Error2_0(channel,snr_ctr,:) = mean(Channel_Est_Error1_0,1);
        Channel_Est_Error2_1(channel,snr_ctr,:) = mean(Channel_Est_Error1_1,1);
        Channel_Est_Error2_2(channel,snr_ctr,:) = mean(Channel_Est_Error1_2,1);
        Channel_Est_Error2_3(channel,snr_ctr,:) = mean(Channel_Est_Error1_3,1);
        Channel_Est_Error2_33(channel,snr_ctr,:) = mean(Channel_Est_Error1_33,1);
        Channel_Est_Error2_4(channel,snr_ctr,:) = mean(Channel_Est_Error1_4,1);
        Channel_Est_Error2_5(channel,snr_ctr,:) = mean(Channel_Est_Error1_5,1);

        
        LA0_2(channel,snr_ctr,:) = mean(LA0_1,1);
        LA1_2(channel,snr_ctr,:) = mean(LA1_1,1);
        LA2_2(channel,snr_ctr,:) = mean(LA2_1,1);
    end

end

Channel_Est_Error3_0 = 10*log10(mean(Channel_Est_Error2_0,1));
Channel_Est_Error3_1 = 10*log10(mean(Channel_Est_Error2_1,1));
Channel_Est_Error3_2 = 10*log10(mean(Channel_Est_Error2_2,1));
Channel_Est_Error3_3 = 10*log10(mean(Channel_Est_Error2_3,1));
Channel_Est_Error3_33 = 10*log10(mean(Channel_Est_Error2_33,1));
Channel_Est_Error3_4 = 10*log10(mean(Channel_Est_Error2_4,1));
Channel_Est_Error3_5 = 10*log10(mean(Channel_Est_Error2_5,1));

LA0_3 = mean(LA0_2,1);
LA1_3 = mean(LA1_2,1);
LA2_3 = mean(LA2_2,1);
%%

for snr_ctr = 1:length(SNR)
    figure;
    subplot(2,1,1)
    plot(reshape(Channel_Est_Error3_0(1,snr_ctr,:),1,L_train-Max_order+1),'b-'); hold all
    plot(reshape(Channel_Est_Error3_1(1,snr_ctr,:),1,L_train-Max_order+1),'k-'); 
    plot(reshape(Channel_Est_Error3_2(1,snr_ctr,:),1,L_train-Max_order+1),'r-');
    plot(reshape(Channel_Est_Error3_3(1,snr_ctr,:),1,L_train-Max_order+1),'b--');
    plot(reshape(Channel_Est_Error3_33(1,snr_ctr,:),1,L_train-Max_order+1),'k--');
    plot(reshape(Channel_Est_Error3_4(1,snr_ctr,:),1,L_train-Max_order+1),'r--');
    plot(reshape(Channel_Est_Error3_5(1,snr_ctr,:),1,L_train-Max_order+1),'g-');
    xlabel('Interation');
    ylabel('MSE')
    legend('Multi-Layered, Bhotto Adaptive-Lambda','Multi-Layered, Paleologu Adaptive-Lambda',...
        'Multi-Layered, Fixed-Lambda','RLS, Bhotto Adaptive-Lambda','RLS, Paleologu Adaptive-Lambda','RLS, Fixed-Lambda','LMS')
%     xlim([4800,6500])

    subplot(3,1,3)
    plot(reshape(LA0_3(1,snr_ctr,:),1,[]),'b-'); hold all
    plot(reshape(LA1_3(1,snr_ctr,:),1,[]),'k-'); 
    plot(reshape(LA2_3(1,snr_ctr,:),1,[]),'r-');
    xlabel('Interation');
    ylabel('L')
    legend('Bhotto Adaptive-Lambda','Paleologu Adaptive-Lambda','Fixed-Lambda')
%     xlim([480,6500])
end

% figure
% plot(abs(e)); hold all
% for c = 1:Max_Layers
%     plot(abs(e_f{c}))
% end












