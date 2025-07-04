%Calculating Bloch_magnetization for MM model 

%assuming a hyperbolic secant inversion pulse profile from Giapitzakis et al MRM 2018
time       =  15*1e-3; %s 
freqOffset = (-2000:10:2000); 
nPoints    = 512;  
timeAxis   =  linspace(0,time,nPoints)*1e+3; 
x = linspace(-pi,pi,nPoints); 
b = 1.904; 
AM = sech(b*x); 
FM = 1017*tanh(-x*1.89); 
PM = (2*pi*cumsum(FM)*time/length(FM)); 
 
%color = ['r','b','m','g','k', 'c', 'y', '-.r', '-.g', '-.b'];
color = colormap(hsv(10));
ih = [0 1 2 -2 -1]; 
for i = 1:length(ih) 
 Bampl      = (1+ih(i)*0.15)*15 *1e-6; %in T 
end
B  = Bampl*AM.*exp(1i*PM); 
 
[alpha, beta] = fSLR( B', freqOffset, time/length(B),0); 
 
%Magnetization vectors 
M0      = zeros(length(freqOffset),3); 
M0(:,3) = 1; 
M = applyRotation( alpha,beta, M0);

Bampl      = 24 *1e-6; %in T 
B          = Bampl*AM.*exp(1i*PM); 
 
TE = 20.1; %ms  %Landheer et al
TR = 2000; %ms  %Landheer et al
[alpha, beta] = fSLR( B', freqOffset, time/length(B),0);
 
%T1 = [280 1030 1513 1746 1777]; % MM Cr-CH2 Cho Cr-CH3 NAA (Values for metabolites are from Deelchand et al at 9.4 T and MM is from Murali-Manohar et al at 9.4 T) 
%T2 = [23.88 81.82 90.11 100.21 110.49]; % (Values are from Murali-Manohar et al at 9.4 T) 

%GM MMs
T1 = [290, 309, 309, 225, 247, 263, 400, 400, 400, 400];
T2 = [27, 39, 18, 17, 14, 20, 26, 21, 18, 21];

%WM MMs
T1 = [284, 280, 270, 253, 319, 313, 383, 518, 379, 434];
T2 = [35, 25, 40, 18, 19, 18, 13, 19, 16, 25];
 
for TI1 = [920]    %For DIR inversion times from Landheer et al      
    for TI2 = [330]   
        for i = 1:length(T1)   
            M = applyRotation( alpha, beta, M0); 
            Mrelax = applyRelaxation( M,T1(i),T2(i),TI1); 
            Mrelax = applyRotation( alpha,beta, Mrelax ); 
            Mrelax = applyRelaxation( Mrelax,T1(i),T2(i),TI2 ); 
            mag.Data{i} = min(Mrelax(:,3)); 
            plot(freqOffset,(Mrelax(:,3)), 'Color', [color(i,:)], 'linewidth',1)
            set(gca, 'fontsize', 12, 'fontweight', 'bold');
            hold on
            set(gca,'xDir','reverse');
            legend('M1', 'M2','M3','M4','M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'fontweight', 'bold')
            xlabel('Frequency (Hz)','fontweight', 'bold', 'fontsize',15);
            ylabel('M_z/M_0', 'fontweight', 'bold', 'fontsize',15);
            axis([-2000 2000 -1 1]); 
        end
    end
end
