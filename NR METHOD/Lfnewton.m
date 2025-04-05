%   Power flow solution by Newton-Raphson method
%   Copyright (c) 1998 by  H. Saadat
ns=0; ng=0; Vm=0; delta=0; yload=0; deltad=0; % Initialize variables to zero  STEP 1.1
% ns: Counter for slack buses, ng: Counter for generator buses,
% Vm: Voltage magnitude array, delta: Voltage angle array,
% yload: Complex load admittance, deltad: Voltage angle in degrees.

nbus = length(busdata(:,1)); % Determine the number of buses from the busdata matrix. STEP 1.2
for k=1:nbus % Loop through each bus to initialize its parameters.
n=busdata(k,1);     % Extract the bus number from the busdata matrix. STEP 2.1
kb(n)=busdata(k,2); Vm(n)=busdata(k,3); delta(n)=busdata(k, 4);  % Assign bus type, voltage magnitude, and angle from busdata. STEP 2.2
Pd(n)=busdata(k,5); Qd(n)=busdata(k,6); Pg(n)=busdata(k,7); Qg(n) = busdata(k,8);     % Extract power demand (Pd, Qd) and power generation (Pg, Qg).
Qmin(n)=busdata(k, 9); Qmax(n)=busdata(k, 10);
Qsh(n)=busdata(k, 11);     % Extract generator reactive power limits (Qmin, Qmax) and shunt admittance (Qsh).
    if Vm(n) <= 0  Vm(n) = 1.0; V(n) = 1 + j*0;       % If voltage magnitude is not specified, assume it to be 1.0 p.u. STEP 2.2
    else delta(n) = pi/180*delta(n); %STEP 2.2
         V(n) = Vm(n)*(cos(delta(n)) + j*sin(delta(n)));    % Convert voltage angle to radians and calculate the complex voltage V(n). STEP2.2
         P(n)=(Pg(n)-Pd(n))/basemva;        % Normalize real and reactive power for per unit system. STEP 2.2
         Q(n)=(Qg(n)-Qd(n)+ Qsh(n))/basemva; % STEP2.2
         S(n) = P(n) + j*Q(n);   % Calculate the complex power S(n) at the bus. STEP 2.2
    end
end
for k=1:nbus  % Loop to count the number of slack and generator buses. STEP 3.1
if kb(k) == 1, ns = ns+1;      % Increment slack bus counter. STEP 3.1
else, end 
if kb(k) == 2 ng = ng+1;   % Increment generator bus counter.  STEP 3.1
else, end
ngs(k) = ng;
nss(k) = ns;   % Store the current count of generator and slack buses for each bus. STEP 3.2
end
Ym=abs(Ybus); t = angle(Ybus);  % Extract the magnitude and angle of the admittance matrix Ybus.
m=2*nbus-ng-2*ns;  % Calculate the total number of unknowns in the system. STEP 4.1
maxerror = 1; converge=1;
iter = 0;  % Initialize error tolerance and iteration counter. STEP 4.2

% Start of iterations
clear A  DC   J  DX  % Clear variables used in the iterative process.
while maxerror >= accuracy & iter <= maxiter % Test for max. power mismatch
for i=1:m
for k=1:m
   A(i,k)=0;      %Initializing Jacobian matrix A with zeros for each iteration STEP 5.1
end, end
iter = iter+1;       % Increment the iteration counter STEP 5.1
for n=1:nbus           % Loop through each bus to calculate the Jacobian matrix and power mismatches. STEP 5.2
nn=n-nss(n);  
lm=nbus+n-ngs(n)-nss(n)-ns;  % Calculate indices for the Jacobian matrix. STEP 5.2
J11=0; J22=0; J33=0; J44=0;   % Initialize elements of the Jacobian matrix for each bus. STEP 5.2 
   for i=1:nbr     % Loop through each branch to calculate the off-diagonal elements of the Jacobian. STEP 5.2
     if nl(i) == n | nr(i) == n   % Check if the branch connects to the current bus. STEP 5.2
        if nl(i) == n,  l = nr(i); end
        if nr(i) == n,  l = nl(i); end                  % Determine the other bus connected by the branch. STEP 5.2
        J11=J11+ Vm(n)*Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
        J33=J33+ Vm(n)*Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));                  % Calculate partial derivatives for the Jacobian matrix elements. STEP 5.2
        if kb(n)~=1
        J22=J22+ Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));
        J44=J44+ Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));                  % Update the Jacobian matrix elements for non-slack buses. STEP 5.2
        else, end
        if kb(n) ~= 1  & kb(l) ~=1    % Only consider non-slack buses for off-diagonal elements. STEP 5.2
        lk = nbus+l-ngs(l)-nss(l)-ns;
        ll = l -nss(l);                       % Determine indices for the off-diagonal elements off diagonalelements of J1  STEP 5.2
        A(nn, ll) =-Vm(n)*Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));    % Update J1 off-diagonal element.  STEP 5.2
              if kb(l) == 0  % update off diagonal elements of J2
              A(nn, lk) =Vm(n)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));end
              if kb(n) == 0  % update off diagonal elements of J3  STEP 5.2
              A(lm, ll) =-Vm(n)*Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n)+delta(l)); end
              if kb(n) == 0 & kb(l) == 0  % update off diagonal elements of  J4  STEP 5.2
              A(lm, lk) =-Vm(n)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));end
        else end
     else , end
   end
   Pk = Vm(n)^2*Ym(n,n)*cos(t(n,n))+J33;
   Qk = -Vm(n)^2*Ym(n,n)*sin(t(n,n))-J11;      % Calculate active and reactive power at the current bus.  STEP 5.2
   if kb(n) == 1 P(n)=Pk; Q(n) = Qk; end   % Swing bus P (% Assign power values for slack bus.)
     if kb(n) == 2  Q(n)=Qk;
         if Qmax(n) ~= 0   % For generator buses, check reactive power limits.   STEP 5.2
           Qgc = Q(n)*basemva + Qd(n) - Qsh(n);
           if iter <= 7                  % Between the 2th & 6th iterations
              if iter > 2                % the Mvar of generator buses are
                if Qgc  < Qmin(n),       % tested. If not within limits Vm(n)
                Vm(n) = Vm(n) + 0.01;    % is changed in steps of 0.01 pu to
                elseif Qgc  > Qmax(n),   % bring the generator Mvar within
                Vm(n) = Vm(n) - 0.01;end % the specified limits.
              else, end
           else,end
         else,end
     end
   if kb(n) ~= 1
     A(nn,nn) = J11;  %diagonal elements of J1 % Update Jacobian matrix and power mismatch vector for active power.  STEP 5.2
     DC(nn) = P(n)-Pk;
   end
   if kb(n) == 0
     A(nn,lm) = 2*Vm(n)*Ym(n,n)*cos(t(n,n))+J22;  %diagonal elements of J2  STEP 5.2
     A(lm,nn)= J33;        %diagonal elements of J3  STEP 5.2
     A(lm,lm) =-2*Vm(n)*Ym(n,n)*sin(t(n,n))-J44;  %diagonal of elements of J4  STEP 5.2
     DC(lm) = Q(n)-Qk;          % Update Jacobian matrix and power mismatch vector for reactive power.  STEP 5.2
   end
end
DX=A\DC';   % Solve the linear system to find the corrections to voltage magnitude and angle.  STEP 5.2
for n=1:nbus  % Apply the corrections to the bus voltages.  STEP 5.2
  nn=n-nss(n);
  lm=nbus+n-ngs(n)-nss(n)-ns;
    if kb(n) ~= 1
    delta(n) = delta(n)+DX(nn); end           % Update voltage angle for non-slack buses. APPLY CORRECTIONS
    if kb(n) == 0
    Vm(n)=Vm(n)+DX(lm); end               % Update voltage magnitude for load buses. APPLY CORRECTION  STEP 5.2
 end
  maxerror=max(abs(DC));     % Check the maximum power mismatch to determine if the solution has converged.  STEP 5.2
     if iter == maxiter & maxerror > accuracy  
   fprintf('\nWARNING: Iterative solution did not converged after ')  %STEP 5.2
   fprintf('%g', iter), fprintf(' iterations.\n\n') % STEP 5.2
   fprintf('Press Enter to terminate the iterations and print the results \n')
   converge = 0; pause, else, end
   
end
 % If the maximum number of iterations is reached without convergence, warn the user.

if converge ~= 1
   tech= ('                      ITERATIVE SOLUTION DID NOT CONVERGE'); else, 
   tech=('                   Power Flow Solution by Newton-Raphson Method');
end 
% Determine the final status of the solution. STEP 6
V = Vm.*cos(delta)+j*Vm.*sin(delta);   % Calculate the final bus voltages.  STEP 6.1
deltad=180/pi*delta;      % Convert voltage angles back to degrees.  STEP 6.1
i=sqrt(-1);     % Define the imaginary unit. STEP 6.1
k=0;  % Initialize counter for slack and generator buses.
for n = 1:nbus    % Loop through each bus to finalize power calculations.
     if kb(n) == 1
     k=k+1;
     S(n)= P(n)+j*Q(n);
     Pg(n) = P(n)*basemva + Pd(n); %STEP 6.2
     Qg(n) = Q(n)*basemva + Qd(n) - Qsh(n); %STEP 6.2
     Pgg(k)=Pg(n); %STEP 6.2
     Qgg(k)=Qg(n);     %june 97   % For slack and generator buses, calculate the total power generation.
     elseif  kb(n) ==2
     k=k+1;
     S(n)=P(n)+j*Q(n);
     Qg(n) = Q(n)*basemva + Qd(n) - Qsh(n);
     Pgg(k)=Pg(n);
     Qgg(k)=Qg(n);  % June 1997 
  end
yload(n) = (Pd(n)- j*Qd(n)+j*Qsh(n))/(basemva*Vm(n)^2);  % Calculate the load admittance at each bus.
end
busdata(:,3)=Vm'; busdata(:,4)=deltad';  % Update the busdata matrix with final voltage magnitudes and angles.  %STEP 6.3
Pgt = sum(Pg);  Qgt = sum(Qg); Pdt = sum(Pd); Qdt = sum(Qd); Qsht = sum(Qsh);
% Calculate the total generated and demanded power. %STEP 7.1

%clear A DC DX  J11 J22 J33 J44 Qk delta lk ll lm   % Clear variables no longer needed.
%clear A DC DX  J11 J22 J33  Qk delta lk ll lm      % Clear variables no longer needed.
