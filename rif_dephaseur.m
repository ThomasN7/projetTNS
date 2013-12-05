function [Brif1,Arif1]=rif_dephaseur(phi1,ob1,W1p)

for k=1:size(W1p,1)
    phirif1(k,1)=k*(min(phi1)-max(phi1))/(size(W1p,1))+max(phi1)-phi1(k,1);
end
    
Hrif1=exp(sqrt(-1)*phirif1);
[Brif1,Arif1]=invfreqz(Hrif1,W1p,ob1,0);