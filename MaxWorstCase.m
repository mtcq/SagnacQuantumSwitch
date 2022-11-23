%Script finds the maximal worst case winning probability for Sequential
%strategies via SDP

%Requires: DefineSets_com_anticom.m, from mtcq

%Author: Marco Túlio Quintino, https://github.com/mtcq, mtcq.mm@gmail.com
%Last update: 11/23/2022

clear all

%Set basic variables
d=2; %dimension
DIM=[d d d d]; %vector with the dimension of spaces Input1 Output1 Input2 Output2

DefineSets_com_anticom;
N=size(Gm,3)+size(Gp,3);

%%%%%%%%% START: dual SDP for the SEQ case %%%%%%%%%
cvx_begin SDP
variable Tp(d^4 d^4) complex semidefinite
variable Tm(d^4 d^4) complex semidefinite
expression pm
expression pp

for im=1:size(Gm,3)
    pm(im)=real(trace(Tm*Gm(:,:,im)));
end
for ip=1:size(Gp,3)
    pp(ip)=real(trace(Tp*Gp(:,:,ip)));
end

Tp+Tm == TR(Tp+Tm,[4],DIM)
TR(Tp+Tm,[3 4],DIM) == TR(Tp+Tm,[2 3 4],DIM)

trace(Tp+Tm)==d^2;

maximise min([pp pm])

cvx_end
%%%%%%%%% END: dual SDP for the SEQ case %%%%%%%%%

p_ij=[pp pm]
MaxWorstCaseSequential=min([pp pm])
