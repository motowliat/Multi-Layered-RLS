clear
M = 10;
e = 0:.001:(2/M);

ro = 1-2*e+(e.^2)*M;
nois = e*M.*ro./(2-e*M)+(1-e*M).^2;

plot(e,ro); hold all
plot(e,nois,'-o')