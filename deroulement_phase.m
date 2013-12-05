function [H1,phi1,W1p]=deroulement_phase(H1,W1)

for k=1:size(H1,1)
    if abs(H1(k,1))<=1e-14
        H1(k,1)=0;
    end
end

F1=find(H1);

phi1=unwrap(angle(H1));

for k=1:(size(phi1,1)-1)
    if phi1(k+1,1)>phi1(k,1)
        phi1(k+1:end,1)=phi1(k+1:end,1)-pi;
    end
end

z1=F1(1,1)-1;
    z2=F1(end,1)+1;

for k=z1:-1:1
    delta=abs(phi1(k+1,1)-phi1(k+2,1));
    phi1(k,1)=phi1(k+1,1)+delta;
end

for k=z2:1:size(phi1,1)
    delta=abs(phi1(k-2,1)-phi1(k-1,1));
    phi1(k,1)=phi1(k-1,1)-delta;
end

W1p=2*pi*W1;
    