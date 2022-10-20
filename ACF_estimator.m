function V = ACF_estimator(X,L_U)

L = L_U;% the length of one slices of X

X = X(1:length(X)-rem(length(X),L_U));% remove extra data
n = length(X)/L;% number of slices

nn = 30; % the number of pieces in ecah slice
LL = floor(L/nn); % length of each piece

for i = 1:n
    y = X((i-1)*L+1:i*L);% slice
    y_middle = y(floor(nn/2)*LL+1:(floor(nn/2)+1)*LL);% the piece in the middle of y
    
    for j = 1:L-LL+1
        v(i,j) = y(j:j+LL-1)'*y_middle;
    end
end
V = mean(v,1);
V = V/max(abs(V));
Add_zero = (L+1)/2-(floor(nn/2)*LL+1);% to have the same peak location with the original ACF
V = [zeros(1,Add_zero) , V(1:L-LL+1-Add_zero)];


