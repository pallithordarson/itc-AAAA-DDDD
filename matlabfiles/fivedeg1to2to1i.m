function b=fivedeg1to2to1i(Ltot,htot,k11,k12,k21);%this is a subfunction within the
%uv1to2b and nmr1to2forb function, it calculates b = concentration of free ligand
%The inputs are Ltot = total conc. of ligand (guest) and htot = total
%conc. of host, K11 and K12, the stepwise binding constants for 1:2 complex

r = size(Ltot,1);%determines the number of datapoints
uu = ones(r,1);%creates a vector full of ones (1) of the same length as data
%this is done to reserve memory and speed up calculation
%Next the coeffiencts for the cubic equation to be solved are determined,
%the cubic equation is a modified from K. A. Connors in Binding Constants, John
%Wiley and Sons, New York, 1987, p. 161 eq. 4.29. See also Tsukube, H.; Furuta, 
%H.; Odani, A.; Takeda, Y.; Kudo, Y.; Liu, Y.; Sakamoto, H.; Kimura, K 
%in Comprehensive Supramolecular Chemistry, ed. Atwood, J. L.; Davis, 
%J. E.D.; MacNicol, D. D.; Vögtle, F.; Lehn, J-M.; Elsevier Science, Oxford, %
%New York, Tokyo, 1996, Vol 8. pp 435.
%The vector uu, as explained above, has the value 1 for all points.
%The equation becomes: a1*x^3 + a2*x^2 + a3*x + a4 = 0 (eq. 1)
a1 = (uu.*(-3.*k11.^2.*k21.^2));
a2 = (uu.*((-4.*k11.^2.*k21)-(6.*Ltot.*k11.^2.*k21.^2)+(3.*htot.*k11.^2.*k21.^2)));
a3 = (uu.*((-k11.^2)+(4.*k11.*k12)-(2.*k11.*k21)-(5.*Ltot.*k11.^2.*k21)+(4.*htot.*k11.^2.*k21)));
a4 = (uu.*((-Ltot.*k11.^2)+(htot.*k11.^2)+(4.*Ltot.*k11.*k12)-(8.*htot.*k11.*k12)+(2.*Ltot.*k11.*k21)+(2.*htot.*k11.*k21)));
a5 = (uu.*(1+(Ltot.*k11)+(Ltot.^2.*k11.*k12)-(4.*Ltot.*htot.*k11.*k12)+(4.*htot.^2.*k11.*k12)));
a6 = (uu.*(-htot));
%a is the matrix were the rows corresponds to the data points and
%the column to the coeffients in the cubic equation
a = [a1 a2 a3 a4 a5 a6];
r = size(a,1);%determs the size of a
b = zeros(r,1);%creates a vector the size of r full of zeros (to reserve memory)
b(1) = htot(1);%sets the first solution (at Ltot=0) =1e-12 since b = 0 could 
%cause problems due to x/0 errors.
for n = 2:r;%starts a loop which solves the cubic equation row by row
   x = roots(a(n,:));%for each row, x = the three solution of eq. 1
   tt=real(x(5));
   b(n) = tt;
end
