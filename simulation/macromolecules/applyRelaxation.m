function M = applyRelaxation( M0,T1,T2,time )  
M           =  zeros(length(M0),3); 
M(:,1)      =   M0(:,1)*exp(-time/T2); 
M(:,2)      =   M0(:,2)*exp(-time/T2); 
M(:,3)      =  (M0(:,3)-1)*exp(-time/T1)+1; 
end 
