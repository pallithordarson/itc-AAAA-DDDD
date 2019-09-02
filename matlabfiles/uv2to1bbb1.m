function b=uv2to1bbb1(Ltot,htot,K11,K12);%this is a subfunction within the
%uv1to2b and nmr1to2b function, it calculates b = concentration of free ligand
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
%New York, Tokyo, 1996, Vol 8. pp 435
%The vector uu, as explained above, has the value 1 for all points.
%The equation becomes: a1*x^3 + a2*x^2 + a3*x + a4 = 0 (eq. 1)
a1 = (uu.*(K11.*K12));
a2 = (uu.*((2.*K11.*K12.*htot)+K11-(Ltot.*K11.*K12)));
a3 = (uu.*((1+(K11.*htot)-(K11.*Ltot))));
a4 = (uu.*(-1.*Ltot));
%a is the matrix were the rows corresponds to the data points and
%the column to the coeffients in the cubic equation
a = [a1 a2 a3 a4];
r = size(a,1);%determs the size of a
b = zeros(r,1);%creates a vector the size of r full of zeros (to reserve memory)
%b(1) = 1e-12;%sets the first solution (at Ltot=0) =1e-12 since b = 0 could 
%cause problems due to x/0 errors.
for n = 1:r;%starts a loop which solves the cubic equation row by row
   x = roots(a(n,:));%for each row, x = the three solution of eq. 1
   if size(x)~=3
       x5 = repmat(x(1)*2,1,3);
       x5(1)=x(1);
       x=x5';
   end
      if isreal(x(1,:));%next, the program looks for the smallest real solution
         w = [x(1,:)];%to eq. 1. Note that some of the solution can be 
      else;%complex numbers, they are first rejected by checking if they are
         w = [0]; %real numbers (the isreal function). 
      end;%if the solutions are not real, they are assign the value of 0
      if isreal(x(2,:));%the three different solution to eq. 1 are given the
         q = [x(2,:)];%values w,q and o after complex solutions have been 
      else;%replaced by the value 0
         q = [0];
      end;
      if isreal(x(3,:));
         o = [x(3,:)];
      else;
         o = [0];
      end;
      t = [w;q;o];%t is the "new" solution vector after complex numbers have been
      if all(t <= 0);%assign the value of 0. Next the program checks if any of the
         b(n) = Ltot(n,:).*rand(1);%solutions are equal to or greater than 0, if 
      else;% not, the value Ltot(n)*random number between 0 and 1 is assigned
         %to that point in the solution vector b. But if one or more solutions 
         %are equal to or greater than 0, the program picks the
         %b(n) = min(t(find(x>0)));%smallest solution which is bigger than 0,   
         b(n) = min(t(x>0));%smallest solution which is bigger than 0,          
         hh = isempty(b(n)); %however they are all = 0, b(n)
         if hh;% (n = n-th data point row) will be an
            %empty vector, this gives hh = TRUE and will call again for a random
            b(n) = Ltot(n,:).*rand(1);%solution of Ltot(n)*rand(1) = random 
         end;%number between 0 and 1.
      end;
      if b(n) > 1111*Ltot(n);%this part of the loops makes sure that if for some
         b(n) = ((Ltot(n))-(Ltot(n))*((K11./(1+K11))));%reason a solution "slips"
      end% through which is much (1111) bigger than Ltot(n), it will be reduced 
   end;%to the value of Ltot(n) - Ltot(n)*a factor given by K11/(1+K11) (a very
   %small number for big K11)
for n = 2:r;%this loop makes a similar check for those datapoints which have 
  if Ltot(n) < 2.*htot(n)%lower total ligand conc (Ltot). than 2*total host
     if b(n) > 1.1*Ltot(n);%conc (htot). First it is checked if the calc. free
        %ligand conc (b(n) is greater than 110% of the total ligand conc (Ltot)
        %If so, the value of b(n) is replaced by half the value of the 
        %total ligand conc.(Ltot);
       b(n) = ((Ltot(n)*0.5));
      end
  else
      if b(n) > 1.1*Ltot(n);%else if the total ligand conc. is greater than
         %2 * total host conc. (htot), the value of b(n) is replaced by
         b(n) = ((Ltot(n)+htot(n))*0.5);%Total ligand conc(Ltot) - half the total
         %host conc. (htot).
     end%Note, these "check" loops are necessary in the beginning of an 
  end%optiminization since incorrect guesses can create quite
  %strange solutions to eq. 1. As the parameter guesses become better, the, it 
end% is more likely that one of the solutions to eq. 1 is resonable.
%ones the end of these loops have been executed, a value b(n) is returned by
%the function uv1to2bb to the "mother" function uv1to2b
