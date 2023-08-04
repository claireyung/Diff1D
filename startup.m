%---------------------------------------------------------------------
%---------------------------------------------------------------------
%
% Startup script for plot settings in MATLAB 
%
%---------------------------------------------------------------------
%
% Dependencies; none. 
%
% Ryans ROMS Matlab and netcdf Utilities 22/2/16
%
%---------------------------------------------------------------------
%---------------------------------------------------------------------

%%addpath(genpath(<path to matlab-utilities>))

%Add to path without using sudo:
% $$$ startdir = pwd;
% $$$ tmp = char(userpath);
% $$$ pathdir = tmp(1:(end-1));
% $$$ cd(pathdir)
% $$$ path(pathdef)
% $$$ cd(startdir)

set(0,'DefaultAxesColorOrder',[0 0 0], ...
      'DefaultAxesLineStyleOrder','-|--|:|-.')
set(0,'defaulttextfontsize',15);
set(0,'defaultaxesfontsize',15);
set(0,'DefaultFigureColor',[1 1 1]);
set(0,'defaultFigureRenderer','painters');
set(0,'DefaultTextInterpreter', 'latex')
set(0,'defaultlinelinewidth',1)
set(0,'defaultaxeslinewidth',1)
set(0,'defaultpatchlinewidth',1)
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultColorbarTickLabelInterpreter','latex');