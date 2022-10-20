clc
clear 
% close all


Cs = 2;% constellation size of mubulation (BPSK)

% double rate
M1 = 50;% must be odd. order of h1, the SI channel
M1_fix = 0;% number of paths of h1 which are fixed and time-invariying 

Max_order = M1;% max of channel orders


lambda_max = .99; % max forgetting factor for adaptive lambda
lambda = .99; %forgetting factor for fixed lambda
delta = 0.001;% initialer of RLS

Max_Layers = 5;% for multi-layered


Input_mode = 3 ;%(1) BPSK, (2) AWGN, (3) AR Process


B = 5e3;% Bandwidth 
Change_interval = 5;%the channel is changing after this many samples in single rate 

Energy_h1 = 0;% in (dB)


L_train = 3000;% symbol length for trainung
CHANNEL = 20; % number of channels to repeat the process
SNR = 30;% evaluating RT-SNR in (dB)
ITR = 5; % number of repeating the process for a fixed snr and channel
COT_range = [20:10:100]*1e-3; % coherence time range


pdp_h1 = exp(-.4*(0:M1-1));% pdp of SI channel in double rate (w)  

E1 = 10^(Energy_h1/10); %in (w) 
pdp_h1 = E1/sum(pdp_h1) * pdp_h1;% adjusted pdp



Channel_Est_Error3_1 = zeros(1,length(COT_range));
Channel_Est_Error3_2 = zeros(1,length(COT_range));
Channel_Est_Error3_3 = zeros(1,length(COT_range));
Channel_Est_Error3_4 = zeros(1,length(COT_range));
Channel_Est_Error3_5 = zeros(1,length(COT_range));

    
COT_ctr = 0;
for COT1 = COT_range % coherence time for non-direct paths of the SI channel
    COT_ctr = COT_ctr+1;
    
    U1 = AFC_generator(COT1,B,Change_interval);% autocorrelation function of non-direct paths, SI channel

    % %-- loading the channel gererator matrices
    for channel = 1:CHANNEL
        WW1_change{channel} = 1/sqrt(2)*(randn(M1,length(U1)+floor((L_train)/Change_interval))+1i*randn(M1,length(U1)+floor((L_train)/Change_interval)));
    end
    % save('WW1_change','WW1_change')
    % load('WW1_change','WW1_change')


    Channel_Est_Error2_1 = zeros(CHANNEL,length(SNR));
    Channel_Est_Error2_2 = zeros(CHANNEL,length(SNR));
    Channel_Est_Error2_3 = zeros(CHANNEL,length(SNR));
    Channel_Est_Error2_4 = zeros(CHANNEL,length(SNR));
    Channel_Est_Error2_5 = zeros(CHANNEL,length(SNR));
  
    for channel = 1:CHANNEL



        %%%%%%%%%%%%%% CHANNEL PREPARATION %%%%%%%%%%%%%%%
        W1 = WW1_change{channel};

        %%%% Training channels generation %%%

        W1_int = 1:length(U1); % interval time of W1 that counts

        hH1_train = zeros(L_train,M1);% h1 channel matrix for training
        for l = 1:L_train-1
            if l == 1 || rem(l,Change_interval) == 0 % because l is the double rate counter
               h1_variable = channel_gen(W1(1:M1,W1_int),pdp_h1,U1,M1_fix);
               W1_int = W1_int+1;
            end
            hH1_train(l,:) = h1_variable; % the channel varies just for signle rate not double rate
        end


     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         





     %%           


        snr_ctr = 0;
        for snr = SNR
            snr_ctr = snr_ctr+1;
            noise_power = Energy_h1-snr;% in (dB)
            N0 = 10^(noise_power/10); % in (w)

            Channel_Est_Error1_1 = zeros(ITR,1);
            Channel_Est_Error1_2 = zeros(ITR,1);
            Channel_Est_Error1_3 = zeros(ITR,1);
            Channel_Est_Error1_4 = zeros(ITR,1);
            Channel_Est_Error1_5 = zeros(ITR,1);
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
                [H1_hat_1,e_1,e_f_1,Lambda_1,LAYERS_1] = AdaptiveLambda_MultiLayered3(x1,y,N0,M1,delta,lambda_max,Max_Layers,Input_mode);
                [H1_hat_2,e_2,e_f_2,Lambda_2,LAYERS_2] = FixedLambda_MultiLayered3(x1,y,N0,M1,delta,lambda,Max_Layers,Input_mode);
                [H1_hat_3,e_3,Lambda_3] = AdaptiveLambda_RLS(x1,y,N0,M1,delta,lambda_max); % simple RLS
                [H1_hat_4,e_4,Lambda_4] = FixedLambda_RLS(x1,y,N0,M1,delta,lambda); % simple RLS
                [H1_hat_5,e_5] = LMS(x1,y,M1,.009);% LMS

                Channel_Est_Error1_1(itr) = sum(mean((abs(hH1_train(1000:L_train-Max_order+1,:)-H1_hat_1(1000:L_train-Max_order+1,:))).^2,1)) / E1;
                Channel_Est_Error1_2(itr) = sum(mean((abs(hH1_train(1000:L_train-Max_order+1,:)-H1_hat_2(1000:L_train-Max_order+1,:))).^2,1)) / E1;
                Channel_Est_Error1_3(itr) = sum(mean((abs(hH1_train(1000:L_train-Max_order+1,:)-H1_hat_3(1000:L_train-Max_order+1,:))).^2,1)) / E1;
                Channel_Est_Error1_4(itr) = sum(mean((abs(hH1_train(1000:L_train-Max_order+1,:)-H1_hat_4(1000:L_train-Max_order+1,:))).^2,1)) / E1;
                Channel_Est_Error1_5(itr) = sum(mean((abs(hH1_train(1000:L_train-Max_order+1,:)-H1_hat_5(1000:L_train-Max_order+1,:))).^2,1)) / E1;

                LA1_1(itr) = mean(LAYERS_1);
                LA2_1(itr) = mean(LAYERS_2);

                clc
                fprintf('itr=%g , snr=%g , channel=%g\n',itr,snr,channel);


            end

            Channel_Est_Error2_1(channel,snr_ctr) = mean(Channel_Est_Error1_1,1);
            Channel_Est_Error2_2(channel,snr_ctr) = mean(Channel_Est_Error1_2,1);
            Channel_Est_Error2_3(channel,snr_ctr) = mean(Channel_Est_Error1_3,1);
            Channel_Est_Error2_4(channel,snr_ctr) = mean(Channel_Est_Error1_4,1);
            Channel_Est_Error2_5(channel,snr_ctr) = mean(Channel_Est_Error1_5,1);

            LA1_2(channel,snr_ctr) = mean(LA1_1);
            LA2_2(channel,snr_ctr) = mean(LA2_1);
        end

    end

    Channel_Est_Error3_1(COT_ctr) = 10*log10(mean(Channel_Est_Error2_1,1));
    Channel_Est_Error3_2(COT_ctr) = 10*log10(mean(Channel_Est_Error2_2,1));
    Channel_Est_Error3_3(COT_ctr) = 10*log10(mean(Channel_Est_Error2_3,1));
    Channel_Est_Error3_4(COT_ctr) = 10*log10(mean(Channel_Est_Error2_4,1));
    Channel_Est_Error3_5(COT_ctr) = 10*log10(mean(Channel_Est_Error2_5,1));

    LA1_3(COT_ctr) = mean(LA1_2,1);
    LA2_3(COT_ctr) = mean(LA2_2,1);

end
%%

figure;
subplot(2,1,1)
plot(COT_range,Channel_Est_Error3_1); hold all
plot(COT_range,Channel_Est_Error3_2);
plot(COT_range,Channel_Est_Error3_3);
plot(COT_range,Channel_Est_Error3_4);
plot(COT_range,Channel_Est_Error3_5);
xlabel('COT (s)');
ylabel('MSE')
legend('Multi-Layered, Adaptive-Lambda','Multi-Layered, Fixed-Lambda','RLS, Adaptive-Lambda','RLS, Fixed-Lambda','LMS')
%     xlim([4800,6500])


subplot(4,1,3)
plot(COT_range,LA1_3); hold all
xlabel('COT (s)');
ylabel('L')
legend('Multi-Layered, Adaptive-Lambda')
%     xlim([480,6500])

subplot(4,1,4)
plot(COT_range,LA2_3); hold all
xlabel('Interation');
ylabel('L')
legend('Multi-Layered, Fixed-Lambda')
%     xlim([480,6500])




% figure
% plot(abs(e)); hold all
% for c = 1:Max_Layers
%     plot(abs(e_f{c}))
% end












