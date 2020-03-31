function datename(xpos,ypos)
% datename(xpos,ypos)
% plot a date and name stamp on current graph
% to be used mainly in fig*.m printing m-files
% xpos, ypos are in units of current plot
% (0,0) = bottom lefthand corner
tmp  = sprintf('%s%s%s','[',date,']');
dn(1) = {tmp};
dn(2) = {'[joshua weitz]'};
tmph = text(xpos,ypos,dn,'fontsize',10);
