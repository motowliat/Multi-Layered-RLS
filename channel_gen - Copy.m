%this function generates the channel coeff
function h = channel_gen(W,pdp,U,M_fix)

h = diag(sqrt(pdp))*W*U.';

h(1:M_fix) = sqrt(pdp(1:M_fix));% the first paths that are fixed

h = h.';
    
    


