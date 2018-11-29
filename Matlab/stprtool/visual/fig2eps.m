function fig2eps(fig)
% FIG2EPS Exports figure to EPS file.
% 
% Synopsis:
%  fig2eps(fig)
%
% Description:
%  This function invokes the standard save dialog to ask
%  the user about the name to which the given figure should
%  be exported in the Encapsulated PostScript (EPS) format.
%
% Input:
%  fig [1x1] Handle of a figure.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 20-Feb-2003, VF

[filename, pathname] = uiputfile('*.eps', 'Save Figure to EPS');

if filename~=0,    

   file=strcat(pathname,filename);
   eval(sprintf('print -depsc -noui -f%f %s',fig,file));
end

return;
