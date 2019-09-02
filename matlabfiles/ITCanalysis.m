%% Prof. Pall Thordarson at UNSW, Sydney, Australia - p.thordarson@unsw.edu.au

% ITC data analysis program

% This program is based on a matlab code from Sahu, D.; Bastidas M.; Lawrence,...
% C. W.; Noid, W. G.; Showalter, S. A. Assessing Coupled Protein Folding and 
%...Binding Through Temperature-Dependent Isothermal Titration Calorimetry, 
%...Methods En-zymol., 2016, 567, 23-45.
%
% With the code from Sahu and co-workers being availabile on Github:
% https://github.com/debsahu/ITC-Fitting-Macros/tree/master

% This program needs to have the *.DAT files (from NITPIC or the MicroCal
% Origin model in the working directory. 
% This program reads in the DAT files to crea an 3-dimensional data array
% called "alldata" where 

% The array alldata will contain:

% first dimension = m-rows = data from each titration point (i)

% second  dimension = 6 columns: 1=experimental heat change (DH) based on
% the integration of the thermograph for each injection point i, 2 =
% injection volume of titrant in microL (INJV) for this injection point i
% 3 = concentration of the injected solute = guest in the cell prior to 
% this injection(Xt), that is the concentration at injection i-1. 
% 4 = concentration of "macromolecule" = host in the cell prior
% to this injection (Mt), that is at injection i-1. 5 = Molar ratio of 
% guest to host (Xt/Mt) after injection i (XMt) and 
%6 = normalised heat change in for injection i in calories per mole of 
% injectant added (this is the data that will be fitted).

% third dimension corresponds to n-repeat experiments (here, usually = 3).

%%
% Program first clear the workspace

clear all

% Labels that will incorporated in the output - need to be manually changed
indicator='FigS17';
label='ITC titration data for the titration of 2Na-2a (host, Run 1 = 0.228 mM, Run 2 = 0.0912 mM and Run 3 = 0.11 mM) with 1a-2I (guest, Rum 1 = 2.09 mM, Run 2 = 0.836 mM and Run 3 = 1.14 mM) in DMSO ';

% 
numIters=3; % number of repeat experiments

%% This reads in the DAT files and converts to matlab array alldata

tables = cell(numIters,1);
FileID = dir('*.dat');
opts = detectImportOptions(FileID(1).name);
% This loop reads in each of the DAT files 
for j=1:length(FileID);
Ts = readtable(FileID(j).name,opts,'ReadVariableNames',false);
tables{j}=[Ts];
end;

% This loop puts all the key data into the array alldata
alldata=zeros(30,6,3);
for k=1:length(FileID);
alldata(:,:,k)=table2array(tables{k});
end
clearvars -except alldata FileID indicator label;

%% This next part of the program then goes through the four binding models
% to be analysied, i.e, 1:1, 1:2, 2:1 and 2:1<->1:1<->1:2 
% In each case a different matlab program file with the prefix 

%Runfitbinding- and suffix itc.m is executed.
%Those files in turn will execute a "fminsearch" (Nelder-Mead) optimization
% of the raw data (from column 6, 2nd dimension in alldata) to the model
% in question and the put all the key outputs into a matlab
% structure "Fitallstuff". This structure and a few other indicators and
% the alldata array are then combined and stored in a structure
% with a name that has the prefix "Fitnow" and suffix = XtoY indicating
% which binding model was used. A plot of the data vs fit is also 
% generated and saved as a tiff file. 
% finally the data structure with the "Fitnow" prefix is saved in a 
% matlab "mat" data file with the name that has a prefix "XtoX" (= binding
% model, a middle part that is a variable linked to the raw DAT file
% name an and suffix "xx.mat".  

%% This loop does the 1:1 binding analysis

for i=1:length(FileID);
filo=FileID(i).name;
imp1=alldata(:,:,i);
Runfitbinding1to1itc; % calls matlab file to perform 1:1 fit to data

%Creates a structure "Fitnow1to1 with all the key outputs
Fitnow1to1{i}={i,Fitallstuff,alldata(:,:,i),filo,indicator,label};
printpre=filo(1:end-4);
suffo=strcat('n2t1to1',filo(1:end-4),'x.tif');
print('-dtiff',suffo,'-r600');
clearvars -except alldata FileID Fitnow1to1 indicator label;
end
filo1=FileID(1).name;
name1to1=strcat('1to1',filo1(1:end-5),'xx');
%Saves all the key outputs to a "mat" file with name that has a 
%prefix 1to1, variable middle (corresponding to raw DAT file) and suffix
%"xx.mat"
save(name1to1);


%% This loop does the 1:2 binding analysis

for i=1:length(FileID);
filo=FileID(i).name;
imp1=alldata(:,:,i);
Runfitbinding1to2itc; % calls matlab file to perform 1:2 fit to data

%Creates a structure "Fitnow1to2 with all the key outputs
Fitnow1to2{i}={i,Fitallstuff,alldata(:,:,i),filo,indicator,label};
printpre=filo(1:end-4);
suffo=strcat('n2t1to2',filo(1:end-4),'x.tif');
print('-dtiff',suffo,'-r600');
clearvars -except alldata FileID Fitnow1to2 indicator label;
end
filo1=FileID(1).name;
name1to2=strcat('1to2',filo1(1:end-5),'xx');
%Saves all the key outputs to a "mat" file with name that has a 
%prefix 1to2, variable middle (corresponding to raw DAT file) and suffix
%"xx.mat"
save(name1to2);

%% This loop does the 2:1 binding analysis

for i=1:length(FileID);
filo=FileID(i).name;
imp1=alldata(:,:,i);
Runfitbinding2to1itc;  %calls matlab file to perform 2:1 fit to data

%Creates a structure "Fitnow2to1 with all the key outputs
Fitnow2to1{i}={i,Fitallstuff,alldata(:,:,i),filo,indicator,label};
printpre=filo(1:end-4);
suffo=strcat('n2t2to1',filo(1:end-4),'x.tif');
print('-dtiff',suffo,'-r600');
clearvars -except alldata FileID Fitnow2to1 indicator label;
end
filo1=FileID(1).name;
name2to1=strcat('2to1',filo1(1:end-5),'xx');
%Saves all the key outputs to a "mat" file with name that has a 
%prefix 2to1, variable middle (corresponding to raw DAT file) and suffix
%"xx.mat"
save(name2to1);
 
%% This loop does the 2:1<->1:1<->1:2 binding analysis
for i=1:length(FileID);
filo=FileID(i).name;
imp1=alldata(:,:,i);
Runfitbinding1to2to1itc; %calls matlab file to perform 2:1<->1:1<->1:2 fit to data

%Creates a structure "Fitnow1to2to1 with all the key outputs
Fitnow1to2to1{i}={i,Fitallstuff,alldata(:,:,i),filo,indicator,label};
printpre=filo(1:end-4);
suffo=strcat('n2t1to2to1',filo(1:end-4),'x.tif');
print('-dtiff',suffo,'-r600');
clearvars -except alldata FileID Fitnow1to2to1 indicator label;
end

filo1=FileID(1).name;
name1to2to1=strcat('1to2to1',filo1(1:end-5),'xx');
%Saves all the key outputs to a "mat" file with name that has a 
%prefix 1to2to1, variable middle (corresponding to raw DAT file) and suffix
%"xx.mat"
save(name1to2to1);