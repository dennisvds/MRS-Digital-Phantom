function M = applyRotation(a ,b, M0) 
 
M       =  zeros(length(M0),3); 
M1tmp   =   conj(a).^2 .*M0(:,1)- b.^2.*M0(:,2)+2*conj(a).*b .*M0(:,3); 
M2tmp   =  -conj(b).^2 .*M0(:,1)+ a.^2 .*M0(:,2)+2*conj(b).*a .*M0(:,3); 
M3tmp   =  -conj(a.*b).*M0(:,1)- a.*b.*M0(:,2)+(abs(a).^2-abs(b).^2) .*M0(:,3); 
M(:,1)= M1tmp;  M(:,2)= M2tmp;  M(:,3)= M3tmp; 
end 
