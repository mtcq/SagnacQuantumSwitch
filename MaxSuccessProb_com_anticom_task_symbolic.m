%Script finds the verifies the upper and lower bounds in rigorous way

%Requires: MaxSuccessProb_com_anticom_task.m, define_sets_switch_symbolic.m from mtcq

%Author: Marco Túlio Quintino, https://github.com/mtcq, mtcq.mm@gmail.com
%Last update: 19/08/2022

MaxSuccessProb_com_anticom_task;

disp('We now use the computer assisted proof method rigorously certify the following upper bound for causal strategies')
UpperSEQsym=sym(904)/sym(1000)
UpperSEQfloat=double(UpperSEQsym)

DefineSets_com_anticom_symbolic;


CSEQsym=sym(CSEQ);
CSEQsym=ProjSeqChannel(CSEQsym,DIM,order);
CSEQsym=MakePSD(CSEQsym);
CSEQsym=CSEQsym/(trace(CSEQsym))*d^2;

N=size(Gm,3)+size(Gp,3);

if IsPSDSym(UpperSEQsym*CSEQsym-Gpsum/N) && IsPSDSym(UpperSEQsym*CSEQsym-Gpsum/N)
    disp('The upper bound was certified by a rigorous computer assisted proof method')
else
    disp('The upper bound could not be certified')
end

