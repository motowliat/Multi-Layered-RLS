clear
clc

M = 50;
lambda = 1-.5/M; %forgetting factor
B = 5e3;% Bandwidth 
COT1 = 40e-3;% coherence time 
Change_interval = 1;

Eps = 1-lambda;
ro = 1-2*Eps+(Eps^2)*M;

[y,U] = AFC_generator(COT1,B,Change_interval);% the target ACF and the filter to this end
L = length(y);
x = 1:L;


figure; hold all
plot(x,y)


q1 = round(interp1(y(1:(L+1)/2),x(1:(L+1)/2),.5));
q2 = round(interp1(y((L+1)/2+1:L),x((L+1)/2+1:L),.5));
q = [q1,q2];
% plot(q,y(q),'ro')
NN(1) = floor((q(2)-q(1))/2); % coherence length
N = NN(1);

for l = 2:6
    
    y_left = circshift(y,-N);
%     plot(x,y_left);

    y_right =  circshift(y,N);
%     plot(x,y_right);


    z = 2*y-y_left-y_right;
%     plot(x,z)

       
    Omega = (lambda.^2/ro).^(abs(x-(L+1)/2));
    Omega(1:q1) = (lambda.^2/ro).^(abs(N));
    Omega(q2:L) = (lambda.^2/ro).^(abs(N));
%     plot(x,Omega)
    
    
    y = z.*Omega;
    plot(x,y)
    
    q1 = round(interp1(y(q1:(L+1)/2),x(q1:(L+1)/2),.5));
    q2 = round(interp1(y((L+1)/2+1:q2),x((L+1)/2+1:q2),.5));
    q = [q1,q2];
%     plot(q,y(q),'ro')
    N = floor((q(2)-q(1))/2); % coherence length
    
    NN(l) = N;

end
plot(x,0.5*ones(1,L),'r-.')

ylim([0 1.1])
legend('${\varphi _{(1)}}[m]$','${\varphi _{(2)}}[m]$',...
    '${\varphi _{(3)}}[m]$','${\varphi _{(4)}}[m]$',....
    '${\varphi _{(5)}}[m]$','${\varphi _{(6)}}[m]$','Interpreter','latex')


%%%%%%%%%%%%%%%%%%%%%
f = (1-Eps).^2/ro;
g = log(1/f);
N = NN(1);
NN_theory(1) = N;
for l = 2:6
    N = ceil((N*log(2) / (3*log(2)+g*N)));
    NN_theory(l) = N;
end

[NN ; NN_theory]



%%%%%%%%%%%%%%%%%%%
% The (l+1)th COT based on the origirnal COT
l = 1:3;
Nl1 = 2*(0.7).^(l)*NN(1) ./ (0.7.^(l-1).*(3.^(l)-1)*g*NN(1)+2*2.08.^l);

NN_new = [NN(1) , Nl1]




