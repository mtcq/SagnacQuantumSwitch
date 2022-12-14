%Script finds the maximal success probabilities for Parallel, Sequential, and general strategies 
%with standard SDP methods

%Requires: DefineSets_com_anticom.m, ProjParChannel.m, ProjSeqChannel.m, ProjNSChannel.m from mtcq

%Author: Marco Túlio Quintino, https://github.com/mtcq, mtcq.mm@gmail.com
%Last update: 19/08/2022

clear all

%Set basic variables
d=2; %dimension
DIM=[d d d d]; %vector with the dimension of spaces Input1 Output1 Input2 Output2

DefineSets_com_anticom;
N=size(Gm,3)+size(Gp,3);

%%%%%%%%% START: dual SDP for the PAR case %%%%%%%%%
cvx_begin SDP
variable C(d^4 d^4) complex semidefinite

sum(Gp,3)/N <= C;
sum(Gm,3)/N <= C;

%TR(C,[2 4],DIM) == TR(C,[1 2 3 4],DIM)
C==ProjParChannel(C,DIM);

minimise trace(C)/d^2
cvx_end
%%%%%%%%% END: dual SDP for the PAR case %%%%%%%%%

lambdaPAR=trace(C)/d^2;
CPAR=C/trace(C)*d^2;

%%%%%%%%% START: dual SDP for the SEQ case %%%%%%%%%
cvx_begin SDP
variable C(d^4 d^4) complex semidefinite

sum(Gp,3)/N <= C;
sum(Gm,3)/N <= C;

%TR(C,[4],DIM) == TR(C,[3 4],DIM)
%TR(C,[2 3 4],DIM) == TR(C,[1 2 3 4],DIM)

C==ProjSeqChannel(C,DIM);

minimise trace(C)/d^2
cvx_end
%%%%%%%%% END: dual SDP for the SEQ case %%%%%%%%%
lambdaSEQ=trace(C)/d^2;
CSEQ=C/trace(C)*d^2;

%%%%%%%%% START: dual SDP for the SEQ-2before1 case %%%%%%%%%
%This considers a sequential scenario where the second box is used before the first
cvx_begin SDP
variable C(d^4 d^4) complex semidefinite

sum(Gp,3)/N <= C;
sum(Gm,3)/N <= C;

%TR(C,[2],DIM) == TR(C,[1 2],DIM)
%TR(C,[2 3 4],DIM) == TR(C,[1 2 3 4],DIM)
C==ProjSeqChannel(C,DIM,[2 1]);

minimise trace(C)/d^2
cvx_end
%%%%%%%%% END: dual SDP for the SEQ case %%%%%%%%%
lambdaSEQ2=trace(C)/d^2;
CSEQ2=C/trace(C)*d^2;


%%%%%%%%% START: dual SDP for the GEN case %%%%%%%%%
cvx_begin SDP
variable C(d^4 d^4) complex semidefinite

sum(Gp,3)/N <= C;
sum(Gm,3)/N <= C;
%sum(Gm,3)/N <= C;
%C-sum(Gm,3)/N>=0;
%TR(C,[4],DIM) == TR(C,[3 4],DIM)
%TR(C,[2],DIM) == TR(C,[1 2],DIM)
C==ProjNSchannel(C,DIM);
minimise trace(C)/d^2
cvx_end
%%%%%%%%% END: dual SDP for the GEN case %%%%%%%%%

lambdaGEN=trace(C)/d^2;
CGEN=C/trace(C)*d^2;

pSuccessPAR=lambdaPAR
pSuccessSEQ12=lambdaSEQ
pSuccessSEQ21=lambdaSEQ2
pSuccessGEN=lambdaGEN

% In case using gate 2 before gate 1 is better, set CSEQ as CSEQ2
if pSuccessSEQ21>pSuccessSEQ12
    CSEQ=CSEQ2;
end
if pSuccessSEQ21>pSuccessSEQ12
    CSEQ=CSEQ2;
    order=[2 1];
else
    order=[1 2];
end
