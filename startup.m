
fprintf('Adding cvx path !! \n');
cDirectory = pwd;

cd cvx-w64;
cvx_setup('cvx_license.dat');
cd(cDirectory);

cvx_solver('Sedumi');
cvx_expert('true');
cvx_quiet('true');

lsDir = dir(cDirectory);
[stringIndex, okButton] = listdlg('listString','Enter Workspace Directory :',...
    'SelectionMode','Single','ListString',{lsDir.name},'ListSize',[200, 250]);

if okButton
    wDirectory = sprintf('%s%s%s',cDirectory,filesep,lsDir(stringIndex).name);
    cd(wDirectory);
    addpath(genpath(pwd));
else
    fprintf('No working path added !!');
end

clc;
