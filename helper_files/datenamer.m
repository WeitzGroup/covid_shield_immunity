function datenamer(xpos,ypos,rot)
% datenamer(xpos,ypos,rot)
% plot a date and name stamp on current graph
% to be used mainly in fig*.m printing m-files
% xpos, ypos are in units of current plot
% rot specifies the rotation in degrees
% (0,0) = bottom lefthand corner
dn = ['[',date,' joshua weitz]'];
tmph = text(xpos,ypos,dn,'fontsize',10);
set(tmph,'rotation',rot);
