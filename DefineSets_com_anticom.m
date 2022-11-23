%Script defines the commuting and anticommuting sets considered in the discrimination task

%Requires: ChoiKetBra.m from mtcq

%Author: Marco Túlio Quintino, https://github.com/mtcq, mtcq.mm@gmail.com
%Last update: 19/08/2022

Id=eye(2);
X=[0 1;1 0];
Y=[0 -sqrt(-1);sqrt(-1) 0];
Z=[1 0;0 -1];
XpY=(X+Y)/sqrt(2);
XmY=(X-Y)/sqrt(2);
XpZ=(X+Z)/sqrt(2);
XmZ=(X-Z)/sqrt(2);
YpZ=(Y+Z)/sqrt(2);
YmZ=(Y-Z)/sqrt(2);
ZmY=(Z-Y)/sqrt(2);
ZpY=(Z+Y)/sqrt(2);

G(:,:,1)=Id;
G(:,:,2)=X;
G(:,:,3)=Y;
G(:,:,4)=Z;
G(:,:,5)=XpY;
G(:,:,6)=XmY;
G(:,:,7)=XpZ;
G(:,:,8)=XmZ;
G(:,:,9)=YpZ;
G(:,:,10)=YmZ;

commutation_relation=zeros(10,10); %Variable which stores commutation relations +1 for commutation -1 for anti, 0 for nothing
setCOM=0;
setANTICOM=0;
for i=1:10
    for j=1:10
        if norm(G(:,:,i)*G(:,:,j)-G(:,:,j)*G(:,:,i))<0.000001
            commutation_relation(i,j)=1;
            setCOM=setCOM+1;
            Gp(:,:,setCOM)=kron(ChoiKetBra(G(:,:,i)),ChoiKetBra(G(:,:,j)));
        end
        if norm(G(:,:,i)*G(:,:,j)+G(:,:,j)*G(:,:,i))<0.000001
            commutation_relation(i,j)=-1;
            setANTICOM=setANTICOM+1;
            Gm(:,:,setANTICOM)=kron(ChoiKetBra(G(:,:,i)),ChoiKetBra(G(:,:,j)));
        end
        %         Gp(:,:,i,j)=kron(ChoiKetBra(G(:,:,i)),ChoiKetBra(G(:,:,j)));
        %         Gm(:,:,i,j)=kron(ChoiKetBra(G(:,:,i)*X*G(:,:,i)'),ChoiKetBra(G(:,:,i)*Z*G(:,:,i)'));
    end
end