function [ss, HG, qHG, HGG, ndhCA, SyC, RR]=uv1to1itc(para,initial,DA) %ss is the outcome of this function, it
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
%This sub-program calculates the binding isotherms based on the 
%parameters (K1, K2,...) provided, either directly or more imporantly by the 
%fminsearch optimisation program. 

%The inputs for this function are
%"para", "initial" and "DA" are the parameters to be fitted
%(K1, K2...), total host and guest concentrations (initial) and 
%raw experimentla data (DA) as explained in more details in 
%the function that calls this program (note para = start1 in that program).

%The output of this function are 
%ss = sum of squares, EA = vector or matrix of Y(DeltaHG), Y(DeltaHG2),...
%values, HG (and HG2 etc..) are the calculated HostGuest complex
%concentration and RR = residual matrix. See the function that calls this
%program for more details. 

%%
%This section extracts the relevant values from the initial input arguments
%above.

%Extracts total host and guest concentrations:
injv = initial(:,1);
xt = initial(:,2);
mt = initial(:,3);
xmt = initial(:,4);
%NEXTLINEVARIABLE (Marker for section that is variable between different fitting programs)   
kn=para(1); 
k=10.^kn;

%Extracts the first binding constant = K1 (K11)from the input "para" (start1 in program calling this function).
DH=para(2); 
vzero=0.00142;
V0=vzero;
vadd=injv(1)+injv(2)+(vzero.*1e+6);
SyC = 2;
SyC=SyC*1e-3;

%K13=para(3);%Extracts the first binding constant = K1 (K11)from the input "para" (start1 in program calling this function).
%ENDVARIABLE (Marker for end of section that is variable between different fitting programs)
npd = 2;

htoti=initial(npd+1,3)*1e-3;
Ltoti=initial(npd+1,2)*1e-3;


I0=[injv(npd+1)*1e-6 mt(npd+2)*1e-3 xmt(npd)];

injv=injv(npd+1:length(injv)-1);

xt=xt(npd+2:length(xt));

mt=mt(npd+2:length(mt));

xmt=xmt(npd+1:length(xmt)-1);

%ndh=ndh(npd+1:length(ndh)-1);

injv=injv*1e-6;

mt=mt*1e-3;
xt=xt*1e-3;


htot=mt;
Ltot=xt;

%%
%This section calculates the required concentration(s) of host-guest
%complex(es).


%NEXTLINEVARIABLE (Marker for section that is variable between different fitting programs)  
%This is Eq 13 in P.Thordarson, Chem. Soc. Rev., 2011, Vol 40, p 1305-1323 
%(see also correction on p 5922 of Chem. Soc. Rev., 2011.)
HG=(0.5.*((Ltot+htot+(1./k))-...
   sqrt(((Ltot+htot+(1./k)).^2)-(4.*((Ltot.*htot))))));

HGi=(0.5.*((Ltoti+htoti+(1./k))-...
   sqrt(((Ltoti+htoti+(1./k)).^2)-(4.*((Ltoti.*htoti))))));


%Here HG = concentration of host-guest complex and htot, Ltot and k are
%obtained from above.
HG=real(HG); %checks that the solution of the above equation is real. 
%removes any complex number etc.
indo=find(isnan(HG)); %searches for HG values that are not real numbers
HG(indo)=sqrt(htot(indo).*Ltot(indo)); %assigns square root(htot*Ltot) 
%GG=Ltot-HG;
HGG=vzero.*DH.*HG;
HGGi=vzero.*DH.*HGi;


qfithg(1)=HGG(1)+(I0(1)/V0)*0.5*(HGG(1)+HGGi)-HGGi;
Qv1(1)=qfithg(1);
Qv1(2:length(HGG)+1,:)=HGG;
stup=Qv1(2:end);
stup1=Qv1(1:end-1);
qHG=stup+((injv/V0).*((stup+stup1)./2))-stup1;
qHG(1)=qfithg(1);

%to any HG values that are not real.
%ENDVARIABLE (Marker for end of section that is variable between different fitting programs)

%%

 % EA(1,:)=EAhost;
ndhCA=qHG./(injv.*SyC);
 
 RR=DA-ndhCA; %Calculates the residual matrix RR as the difference between


%UA = raw data matrix and HGG*EA - calculated data which is the matrix 
%product (HGG x EA) of hostguest complex matrix (HGG) and the EA = vector of
%Y(DeltaHG)... values.


%CA=HGG*EA; %Calculates the data again as the matrix product HGG x EA so it can 



%be passed back to the original function as one of the outputs of this
%function.
%ENDVARIABLE (Marker for end of section that is variable between different fitting programs)

%%
%This final segment calculates the sum of squares (ss) which is the target
%function for the optimisation function fminsearch to make as small as
%possible (=0)


rr=RR(:); %Converts the residual matrix RR from above to a column vector = rr


%save running b HGG K11 K12 K13 UA EA;

ss=sum(rr.^2); %calculates the sum of squares of rr (residual matrix).
%end of this function.


