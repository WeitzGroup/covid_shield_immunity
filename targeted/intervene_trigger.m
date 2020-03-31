function [value,isterminal,direction] = intervene_trigger(t,y,pars,agepars)
% function [value,isterminal,direction] = intervene_trigger(t,y,pars,agepars)
% 
% Check to see if total cases exceeds a trigger

Itot = 1-sum(y(agepars.S_ids));
value = Itot - pars.Itrigger;
isterminal = 1;
direction = 0;
