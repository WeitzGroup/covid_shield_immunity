function dydt = sir_shield(t,y,pars)
% function dydt = sir_shield(t,y,pars)
% 
% Shield model

% Variables
S=y(1);
I=y(2);
R=y(3);
dydt=zeros(3,1);

% Model
dydt(1) = -pars.beta*S*I/(S+I+(1+pars.alpha)*R);
dydt(2) = pars.beta*S*I/(S+I+(1+pars.alpha)*R)-pars.gamma*I;
dydt(3) = pars.gamma*I;
