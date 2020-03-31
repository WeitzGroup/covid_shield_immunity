function psprint(filename)
% psprint('filename')
% prints current figure with 
% the following flags
% -deps2
% .ps will be appended
filenameps = sprintf('%s.ps',filename);
print(filenameps,'-deps2');
filenamjpg = sprintf('%s.jpg',filename);
print(filenamjpg,'-djpeg80');
% trouble with the color
%print(filenameps,'-deps');
