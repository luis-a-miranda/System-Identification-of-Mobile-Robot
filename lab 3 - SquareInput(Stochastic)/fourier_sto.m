function [PHASE,GAIN,W]=fourier_sto(ir_cut)

%building axis dimension and vector sizes
W=logspace(0,2); 
preal=zeros(length(W),1);
pimag=zeros(length(W),1);
Gjw=zeros(length(W),1);
PHASE=zeros(length(W),1);
GAIN=zeros(length(W),1);

for m=1:length(W)
    for k=1:length(ir_cut)
        preal(m)=ir_cut(k)*cos(k*(W(m)*(1/30)))+preal(m);
        pimag(m)=ir_cut(k)*sin(k*(W(m)*(1/30)))+pimag(m);
    end
    
Gjw(m)=(1/30)*preal(m)-1j*(1/30)*pimag(m);
PHASE(m)=angle(Gjw(m));
GAIN(m)=20*log10(abs(Gjw(m)));
end
