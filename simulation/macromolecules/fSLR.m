function [alpha, beta] = fSLR( B1, freqOffset, deltaT,rephasePoint) 
%FSLR Summary of this function goes here 
% B1 profile in T 
% freqOffset in Hz 
% deltaT     in secs 
gammaHz   = 42.576 *1e+6; %Hz/T  
gammaRads = gammaHz*2*pi; %rad/(s*T) 
freqOffset = freqOffset.'; 
G         = freqOffset./gammaHz; %Tesla 
a = zeros(length(freqOffset),2);    
b = zeros(length(freqOffset),2); 
a(:,1)      = ones(length(freqOffset),1);   
flagRephase = true; 
for i=1:length(B1) 
    if flagRephase  
        if (i == rephasePoint +1) 
        G= -G;   
        flagRephase = false; 
        end 
    end 
    B1eff = B1(i).*ones(size(freqOffset)); 
    phi = -gammaRads * deltaT * sqrt( abs(B1eff).^2 + G.^2); 
    n0 = (gammaRads*deltaT./abs(phi));   
    n   = [real(B1eff).*n0, imag(B1eff).*n0, G.*n0]; 
    av = cos(phi/2) -1i*n(:,3).*sin(phi/2); 
    bv = -1i*(n(:,1)+1i*n(:,2)).*sin(phi/2);
    a(:,2) = av.*a(:,1) - conj(bv).*b(:,1); 
    b(:,2) = bv.*a(:,1) + conj(av).*b(:,1); 
    a(:,1) = a(:,2); b(:,1) = b(:,2); 
end 
alpha = a(:,2);  beta  = b(:,2); 
end
