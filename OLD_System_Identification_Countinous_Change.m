clc
clear 
% close all


Cs = 2;% constellation size of mubulation (BPSK)

% double rate
M1 = 11;%51;% must be odd. order of h1, the SI channel
M1_fix = 0;% number of paths of h1 which are fixed and time-invariying 

Max_order = M1;% max of channel orders


lambda_max = .99;% forgetting factor of RLS
delta = 0.001;% initialer of RLS
miu_trn = .1;% damping factor for training period
miu = .005; % damping factor for the test period


Layers_Adaptivity = 1;
Max_Layers = 5;


B = 5e3;% Bandwidth 
COT1 = 40e-3;% coherence time for non-direct paths of the SI channel
Change_interval = 5;%the channel is changing after this many samples in single rate 

Energy_h1 = 0;% in (dB)


L_train = 10000;% symbol length for trainung
CHANNEL = 10; % number of channels to repeat the process
SNR = [0:5:55];% evaluating RT-SNR in (dB)
ITR = 5; % number of repeating the process for a fixed snr and channel

pdp_h1 = exp(-.4*(0:M1-1)).*[1 1 0 1 1 0 0 1 0 1 0];% pdp of SI channel in double rate (w)  
pdp_h1 = pdp_h1(1:M1);% pdp of SI channel in double rate (w) 


E1 = 10^(Energy_h1/10); %in (w) 
pdp_h1 = E1/sum(pdp_h1) * pdp_h1;% adjusted pdp


U1 = AFC_generator(COT1,B,Change_interval);% autocorrelation function of non-direct paths, SI channel

% %-- loading the channel gererator matrices
% for channel = 1:CHANNEL
%     WW1{channel} = 1/sqrt(2)*(randn(M1,length(U1)+floor((L_train)/Change_interval))+1i*randn(M1,length(U1)+floor((L_train)/Change_interval)));
% end
% save('WW1','WW1')
load('WW1','WW1')


Channel_Est_Error2 = zeros(1,length(SNR));
Channel_Est_Error_tild2 = zeros(1,length(SNR));
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
        hH1_train(l,:) = h1_variable; % the channel varies just for signle rate not double rate
    end
                  

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
     
 
 

 
 %%           
         
    
    snr_ctr = 0;
    for snr = SNR
        snr_ctr = snr_ctr+1;
        noise_power = Energy_h1-snr;% in (dB)
        N0 = 10^(noise_power/10); % in (w)
        
        Channel_Est_Error1 = zeros(1,ITR);
        Channel_Est_Error_tild1 = zeros(1,ITR);
        for itr = 1:ITR

            
        %%%% Training Phase
        
            X1 = qammod(randi([0,Cs-1],1,L_train),Cs);%
            x1 = [zeros(1,M1-1) , X1];% added zeros symbos to check the BER

            
            s = zeros(1,L_train-Max_order+1);
            for l = 1:L_train
                h1_variable = hH1_train(l,:);
                s(l) = h1_variable*x1(l+M1-1:-1:l).';%SI signal
            end

            
            noise = sqrt(N0/2)*(randn(1,L_train-Max_order+1)+1i*randn(1,L_train-Max_order+1)); %ambiant noise
            
            y = s(1:L_train-Max_order+1)+noise(1:L_train-Max_order+1);% Training with noise

            
            %--- using the signals to train the channels
            [H1_hat,e,e_f,Lambda,LAYERS] = Train_MultiLayered(x1,y,N0,M1,delta,lambda_max,Layers_Adaptivity,Max_Layers);
%             [H1_hat,e,e_f,Lambda,LAYERS] = Train_MultiLayered2(x1,y,N0,M1,delta,lambda_max,Layers_Adaptivity,Max_Layers);

                        
            Channel_Est_Error1(itr) = sum(mean((abs(hH1_train(1000:L_train-Max_order+1,:)-H1_hat(1000:L_train-Max_order+1,:))).^2,1)) / E1;
             
            clc
            fprintf('itr=%g , snr=%g , Layers=%g, channel=%g\n',itr,snr,mean(LAYERS),channel);
            Channel_Est_Error2
                
        end

        Channel_Est_Error2(channel,snr_ctr) = mean(Channel_Est_Error1);
    end

end

Channel_Est_Error3 = mean(Channel_Est_Error2,1);


figure
semilogy(SNR,Channel_Est_Error3); hold all
legend(['LayersAdaptivity:',num2str(Layers_Adaptivity), ',  MaxLayers:',num2str(Max_Layers)])
xlabel('SNR');
ylabel('MSE')

figure
plot(abs(e)); hold all
for c = 1:Max_Layers
    plot(abs(e_f{c}))
end












