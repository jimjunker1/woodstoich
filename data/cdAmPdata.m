%% cdAmPdata
% cd to the AmPdata directory

function WD = cdAmPdata
% created 2020/01/12 by Bas Kooijman

%% Syntax
% WD = <../cdAmPdata.m *cdData*>

%% Description
% cd to the AmPdata directory
%
% Output
%
% * WD: current path 

%% Remarks
% Intended use: WD = cdAmPdata; ..code.. cd(WD)

WD = pwd; path = which('cdAmPdata'); 
if ismac
 ind = strfind(path,'/'); 
else
 ind = strfind(path,'\');
end
cd(path(1:ind(end)));
