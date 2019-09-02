    %%
%               FIT BINDING
%
%
%(C) Dr. Pall Thordarson
%School of Chemistry
%UNSW
%AUSTRALIA
%p.thordarson@unsw.edu.au
%
%Please cite: P.Thordarson, Chem. Soc. Rev., 2011, Vol 40, p 1305-1323 
%when using this program.
%
%A program for determining binding constants from titration experiments in
%supramolecular chemistry
%
%This sub-program prepares the raw data and the initial binding constants for the
%fitting program used, followed by some post-calculations after which all the 
%key data is captured in structur files, and saved.
%

% Clear all the data loads the data needed
%clear all;

%Loads the input file
%imp1=importdata('run1ofs19.txt');
%load s19run1;

%This file includes the required input data:

% DA = raw data to be fitted
% htot = Host concentration
% Ltot = Guest concentration
% initial = [htot Ltot] - maxtrix with the host and guest concentration

% Next we set the initial guesses - this can be edited in the file

start1 = [5 6 -5520 11520];
% First number = K11, second K12 and third = K13

oldstart=start1;

%%
%This segment prepares the data to be fitted

npd=2;

initial=imp1(:,2:5);

DAold=imp1(:,6);
DA=DAold(npd+1:length(DAold)-1);

%%
%This segment execute the fminsearch optimisation process. It is set up as loop 
%that will run up to 10 times or until there is no significant changes in
%the sum of residual squares (ss) between the raw and fitted data which
%fminsearch is trying to make = 0;

options = optimset('TolFun',1e-18,'TolX',1e-18); %Sets the tolerance options for fminsearch

ssold=1e60; %initial (artificially high) value for sum of squares. Stops
%the loop from terminating before the first run.
itn=0; %initial counter for loop function
ssdiff=1e-8;%criteria for loop to stop where ssdiff=ss(this cycle) - ss(last cycle)
tic;%Starts a clock to time the process
while itn<10 %start of loop - itn<10 means it will stop after 10 iterations.

    %The next line executes the fminsearch optimisation.
    %INPUTS for fminsearch are:
    
    %'xx#to#fitbind$$$$' = fitting function that contains the program that calculates
    %naming convention for fitting function:
    %xx = uv or nmr or flu...
    %#to# = stochiometry, e.g. 1to1, 1to2, 2to1,...]\
    %$$$ = optional suffix for variants such as "stat" = statistical 1:2
    %model
    
    %the fitted data (isotherms) and then compares the outcome with the raw
    %data to calculate the sum of squares (ss) which fminsearch tries to
    %minimize.
    
    %start1=input of the initial guess(ss) for the parameters to be fitted
    %(K1, K2....)
    %options = the options for fminsearch as set above.
    %initial = a two column vector with initial host and guest
    %concentrations
    %DA = raw data to be fitted
%NEXTLINEVARIABLE (Marker for section that is variable between different fitting programs)   
[results, fval, exitflag, output] = fminsearch('uv2to1itc',start1,options,initial,DA);
%ENDVARIABLE (Marker for end of section that is variable between different fitting programs)

%OUTPUT from fminsearch
%results = optimised parameters = K1, K2....
%fval = value of sum of squares (ss, target = 0)
%exitflag = messages from fminsearch about success or otherwise
%output = string containing information about the performance of fminsearch
%such as no of iterations

%loop now checks the outcome from fminsearch:

criter=(ssold-fval)/ssold; %measures difference between this cycle and the last
iter(itn+1)=output.iterations; %saves the no of iterations from fminsearch as = iter
if abs(criter) <=ssdiff; %if the ss-difference between this and the last cycle is 
    %less than (<) the preset criteria above
    break % the loop terminates.
end
ssold=fval; %if loop didn't terminate, ssold (previous cycle) is now reassigned 
%as ss from this cycle (fval) 
start1=results; %The initial guess for the next cycle is now changed to the 
%best fit of parameters for this cycle.
itn=itn+1; %increments the number of cycles completed (until 10)
%start1=results;
end %end of loop.

fittime=toc; %records the time the optimisation took as fittime
%%
%Post-processing - part 1 - calculating the final spectra 

tic; %starts a clock to measure the time for all the post-processing processes.

para=results; %takes the best fit from loop agove (results) and assigns it as = para

%Fitting function 'xx#to#fitbind$$$$' is now executed once more with the
%best fit from optimisation loop above.
%Input = para = K1, K2...., initial and DA = initial concentration and raw
%data as before

%NEXTLINEVARIABLE (Marker for section that is variable between different fitting programs)   
[ss, HG, H2G, qHG, HGG, ndhCA, SyC, RR]=uv2to1itc(para,initial,DA);

%ENDVARIABLE (Marker for end of section that is variable between different fitting programs)

results(1:2)=10.^results(1:2);
oldstart(1:2)=10.^oldstart(1:2);



%outputs:

%ss = sum of squares from fitting function
%EA = vector or matrix of Y(DeltaHG), Y(DeltaHG2),... values. 
%E.g. for 1:1 binding and: 
%UV; EA = [DeltaEpsilonWavelength1,DeltaEpsilonWavelength2,....]
%NMR; EA = [DeltadeltainppmorHzProton#1, DeltadeltainppmorHzProton#2,...]
%see e.g. the first term on the right of the "=" in eq. 14
%and 15 in P.Thordarson, Chem. Soc. Rev., 2011, Vol 40, p 1305-1323.
%HG = Calculated concentration of the HG species from the fitting process
%RR = The calculated matrix of residuals (Raw data - Fitted data) from the
%fitting process.

DA=DAold; %removes the "dummy" ttb(2) (see above) again from 
%raw data column DA

%DA=reshape(DA,tth+1,ttb(2)); %Regenerates the original matrix of 
%raw data (DA) as a m x n matrix with tth = number of rows and
%ttb(2) = number of columns
%% 
%Post-processing - part 2 - calculates key statistics 

ndata=length(ndhCA); %measures the number of experimental datapoints in raw data
%ndata = number of titration point x number of spectra in global fit
npara=length(results); %measures the binding of binding constants fitted (usually 1 or 2)
%nEA=1; %measures the number of datapoints in global fit fitted 
%(= no of columns in original DA x number of species).
%nEA2=1;%measures the number of spectra in global fit fitted (= no of columns in original DA).
deFree=ndata-npara; %deFree = degree's of freedom 
SE=(sqrt(ss/(deFree))); %Standard devitation of the y-estimate (Eq 52 in 
%P.Thordarson, Chem. Soc. Rev., 2011, Vol 40, p 1305-1323.
chiSQ=(sqrt(ss/(deFree-1))); %Chi-square = Eq 53 in P.Thordarson, Chem. Soc. Rev., 2011, Vol 40, p 1305-1323.

%The next section is to calculate the uncertainty for the linear least square process
%that was used in the fitting function program to obtain EA (spectral shifts) by matrix
%division (linear regression).

%%
DAA2=DA(:);
%CAA2=CA(:);
sigmas=ss./ndata;
maxloglik=-((ndata./2).*(log(2.*pi)))-((ndata./2).*(log(sigmas)))-((ss./(2.*sigmas)));
myBIC =(-2.*maxloglik)+(log(ndata).*(ndata-deFree));
%%

%This is done here with the inverted "curvature" matrix approach as shown
%in Eq 4.38, page 123 in "Practical Data Analysis in Chemistry" by M. Maeder and Y.-M. Neuhold,
%Elsevier, 2007 (Vol 26 in the series Data handling in Science and Technology).
invHG=inv(HG'*HG); %Calculates the inverted curvature matrix from HG = calculated host-guest complex. 
uncertEAinv=sqrt(diag(invHG)); %Takes the square root of the diagonal of the inverted curvature matrix.
uncertEA=SE*uncertEAinv; %calculates the standard deviation of the EA parameters.

%perEA=(uncertEA./EA(:)).*100; %calculates the %standard deviation of each of the EA parameters.

%Next we calculate the covariance of fit (covfit). This parameter is less
%sensitive to no of parameters than chiSQ and SE - see Eq 54 in P.Thordarson, Chem. Soc. Rev., 2011, Vol 40, p 1305-1323.
covRR=cov(RR(:)); %Takes the covariance of the residual matrix after converting it to a column vector
covDA=cov(DA(:));%Takes the covariance of the raw data (-initial values)after converting it to a column vector
covFit=covRR./covDA; %calcules the covfit according to Eq 54 in P.Thordarson, Chem. Soc. Rev., 2011, Vol 40, p 1305-1323.
%

%
%%
%Post-processing - part 3 - prepares parameters for later plotting processes / GUI's

xx=initial(:,2)./initial(:,3); %Calculates the number of Guest added as xx = Total Guest(Ltot) /Total Host (htot) concentration.
dev=RR./SE; %Normalises the residual matrix (RR) by dividing with standard error of the Y-estimate SE. This parameter (dev) is ideal for
%use to do a residual plot (see e.g top of page 1320 of P.Thordarson, Chem. Soc. Rev., 2011, Vol 40, p 1305-1323) 
sdev2=sum(dev,1); %Does the sum of the normalised residuals (dev) column by column (per spectral line in DA)
%Hfree=htot-HG;

plot(xx(npd+1:end-1),ndhCA,'r');
hold on;
plot(xx,DAold,'b*');
hold off
%%
%Post-processing - part 3 - gathers data in structures and saves it
%then savedfor 
timeoffit=datestr(clock, 0); %Creates a date/time stamp for future reference

%Creates structures Fitresults,FitversExpData and FitConcAndSpectra for
%fitting results, Experimental vs. calc. data and calc
%concentrations/spectra, respecively
Fitresults=struct('K1',results,'Initial_guess_K1',oldstart);
%FitversExpData=struct('Equiv_GuestAdd',xx,'ExpData',DA+hostA,'FittedData',CA+hostA,'Residual_Data',RR,'Host_Conc',htot,'Guest_Conc',Ltot,'Wavel_Signals',wavel);

FitConcAndSpectra=struct('Host_Guest_Q_data',[HGG],'calculated_molar_heat',ndhCA,'host_guest_1t1_conc',HG,'heat_per_injection',qHG,'residuals',RR,'SyC',SyC,'host_guest_2t1_conc',H2G);
%The FitStats structure collects various statistical data + the name of the
%program used.
%NEXTLINEVARIABLE (Marker for section that is variable between different fitting programs)   
FitStats=struct('SumofSq_ss',ss,'StandardError_Y_estim_SE',SE,'Chi_squared',chiSQ,'Covariance_Fit_covfit',covFit,...
    'maxlikel',maxloglik,'BIC',myBIC,'sigma',sigmas,'DegreesFreedom',deFree,'Program_Used','fitbinding1to3uv');
%ENDVARIABLE (Marker for end of section that is variable between different fitting programs)

plottime=toc; %stops the clock for measuring the post-processing procedures.

%FittimePerform collects various fitting times, no of iteration data and
%date/time stamp.
FitTimePerform=struct('Time_of_fitting',fittime,'Cycles_Number_iter',[itn output.iterations],'Time_Plotting',plottime,'Date_Time_Stamp',timeoffit,'Other_output',output);


%The following loop creates the final Fitallstuff cell array. 
%The loops add a new line to the preexisting Fitallstuff for each iteration
%cycle as "st" is increased (see home directory program)
st = 1;
Fitallstuff={st,results, Fitresults, FitConcAndSpectra, FitStats, FitTimePerform};
   

%Saves the new Fitallstuff cell array by appending it to the original data file (as
%specified by "str" in home directory program).
%save ('fitted1to3','Fitallstuff','-append');
%end of program.