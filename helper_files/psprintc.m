function psprintc(filename)
% psprint('filename')
% prints current figure with 
% the following flags
% -depsc2
% .ps will be appended
%
% print in color
filenameps = sprintf('%s.pdf',filename);
print(filenameps,'-dpdf');
filenameps = sprintf('%s.ps',filename);
print(filenameps,'-depsc2');
filenamjpg = sprintf('%s.jpg',filename);
print(filenamjpg,'-djpeg80');
