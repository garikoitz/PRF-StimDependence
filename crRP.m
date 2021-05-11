function rootPath=crRP()
%
%        rootPath =mrvRootPath;
%
% Determine path to root of the coverageReading directory
%
% This function MUST reside in the directory at the base of the
% coverageReading directory structure
%
% Wandell

rootPath=which('crRP');

[rootPath,fName,ext]=fileparts(rootPath);

return
