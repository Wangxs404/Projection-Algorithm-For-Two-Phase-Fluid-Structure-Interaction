phiVolEnd=phi ; phiVolEnd(phiVol<0)=0; 
VolumnEnd=sum(sum(phiVol))*(dx*dy);
fprintf('VolumnLeaf= %s %s\n',num2str((VolumnIni-VolumnEnd)/VolumnIni),'%');