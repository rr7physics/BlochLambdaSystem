%Solves Bloch equations for a Lambda system.
%Equations taken from Iyyanki V. Jyotsna and G. S. Agarwal. Phys. Rev. A 52, 3147 (1995)

%Parameters
gamma = 1.;
gammaA = gamma;
gammaB = gamma;
%Detunings
Delta1 = 0;
Delta2 = 0;
%STD of TW Gaussian
sigmaTW = 5.;
%STD of SW Gaussian
sigmaSW = 5.;
%Gaussian peak time for SW
tSW = 20.;
%Gaussian peak time for TW
tTW = 10.;
%Cavity finess
R = 100.;
%Wave-vector assuming lambda=1, so x range is between -0.1 and 0.1
k = 2*pi;

%Amplitude of OmegaTW
OmegaTW0 = linspace(0.25,4.25,100);
%Amplitude of OmegaSW as defined in Mompart et al. (2009)
OmegaSW0 = sqrt(R)*OmegaTW0;

%Space interval being probed
spaceInterval = linspace(-0.1,0.1,1001);
%Evolution time range
tInterval = [0,40];
%Initial populations and coherences, indexing CC, CB, CA, BB, BA, AA
initRho = [0,0,0,0,0,1]';

%Initialization of |a> level population result matrix
rhoAA = [];

%Loops over OmegaTW0 parameters
for i = 1:length(OmegaTW0)
    
    OmegaTW0val = OmegaTW0(i);
    OmegaSW0val = OmegaSW0(i);
    
    if mod(i,5) == 0
        fprintf('Calculated for %d / %d parameters\n',i,length(OmegaTW0))
    end
    
    rhoAAz = [];
    
    %Parallel loop solving Bloch equations for different spatial positions
    parfor j = 1:length(spaceInterval)
      
        z = spaceInterval(j);
        [t,y] = ode45(@(t,rho) blochEquations(t,rho,Delta1,Delta2,gammaA,gammaB, ...
            OmegaTW0val,OmegaSW0val,tTW,tSW,sigmaTW,sigmaSW,k,z) , tInterval , initRho);
        
        rhoAAz = [rhoAAz, y(end,6)];
    end
    
    rhoAA = [rhoAA;rhoAAz];
end

%plot
fig1 = pcolor(spaceInterval,OmegaTW0,rhoAA);
xlabel('x/lambda')
ylabel('OmegaTW0/gamma')
fig1.EdgeColor = 'none';
colorbar;

function dYdt = blochEquations(t,rho,Delta1,Delta2,gamma1,gamma2,OmegaTW0,OmegaSW0,tTW,tSW,sigmaTW,sigmaSW,k,z)
%Index order CC, CB, CA, BB, BA, AA -> 11,12,13,22,23,33 -> 1,2,3,4,5,6

OmegaTW = OmegaTW0*exp(-(t-tTW)^2/sigmaTW^2);
OmegaSW = OmegaSW0*sin(k*z)*exp(-(t-tSW)^2/sigmaSW^2);

dYdt = [ -2*(gamma1+gamma2)*rho(1)+1j*OmegaSW*conj(rho(3))+1j*OmegaTW*conj(rho(2))-1j*conj(OmegaSW)*rho(3)-1j*conj(OmegaTW)*rho(2);
         -(gamma1+gamma2-1j*Delta2)*rho(2)+1j*OmegaSW*conj(rho(5))+1j*OmegaTW*(rho(4)-rho(1));
         -(gamma1+gamma2-1j*Delta1)*rho(3)+1j*OmegaTW*rho(5)+1j*OmegaSW*(1-2*rho(1)-rho(4));
         2*gamma2*rho(1)-1j*OmegaTW*conj(rho(2))+1j*OmegaTW*rho(2);
         1j*(Delta1-Delta2)*rho(5)-1j*OmegaSW*conj(rho(2))+1j*conj(OmegaTW)*rho(3);
         2*gamma1*rho(1)-1j*OmegaSW*conj(rho(3))+1j*conj(OmegaSW)*rho(3)];
end


