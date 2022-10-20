function [V_shift,U] = AFC_generator(COT,B,a)

betta = log(2)./COT;
t = 0:a/B:-log(1e-3)/betta;%step size is a points

V1 = exp(-betta*t);
V2 = V1(end:-1:2);
V = [V1 , V2];% this is the auto correlation

U = fftshift(ifft(sqrt(real(fft(V))))); %filter to biuld the function
U = U/sqrt(sum(U.^2));% normalized energy
% plot(fftshift(V)); hold all
% plot(U)
V_shift = fftshift(V);% target ACF