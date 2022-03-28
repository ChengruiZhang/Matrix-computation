% Julius_Duthoy_IPM_test.m
function Julius_Duthoy_IPM_test
%Converging Example
A=[-4 14 0;-5 13 0;-1 0 2];
x0=[1;1;1];
M=1000;
val=eig(A);
[V,D]=eig(A);
[h,k]=size(V);
for i=1:k
   V(1:h,i)=V(1:h,i)/norm(V(1:h,i),inf); 
end
disp('Actual Eigenvalues of A')
disp(val)
disp('Actual Corresponding Unit Eigenvectors of A')
disp(V)
[ lambda,x ] = InversePM( A,x0,M )


%Non-Converging Example
%A=[2 1 -1;-1 0 2;-1 -2 4]
%A=[7 0 -3;-9 -2 3;18 0 -8];
A=[-9 11 -21 63 -252;70 -69 141 -421 1684;-575 575 -1149 3451 -13801;3891 -3891 7782 -23345 93365;1024 -1024 2048 -6144 24572];
B=inv(A);
x0=[1;1;1;1;1];
val=eig(A);
index=find(val==min(val));
Magnitude_Min_Eigenvalue=abs(val(index))
[V,D]=eig(A);
[h,k]=size(V);
for i=1:k
   V(1:h,i)=V(1:h,i)/norm(V(1:h,i),inf); 
end
disp('Actual Eigenvalues of A')
disp(val)
disp('Actual Unit Eigenvectors of A')
disp(V)
[ lambda,x ] = InversePM( A,x0,M )

%Non-Converging Example
%A=[2 1 -1;-1 0 2;-1 -2 4]
A=[7 0 -3;-9 -2 3;18 0 -8];
%A=[-9 11 -21 63 -252;70 -69 141 -421 1684;-575 575 -1149 3451 -13801;3891 -3891 7782 -23345 93365;1024 -1024 2048 -6144 24572];
B=inv(A);
x0=[1;1;1];
val=eig(A);
index=find(val==min(val));
Magnitude_Min_Eigenvalue=abs(val(index))
[V,D]=eig(A);
[h,k]=size(V);
for i=1:k
   V(1:h,i)=V(1:h,i)/norm(V(1:h,i),inf); 
end
disp('Actual Eigenvalues of A')
disp(val)
disp('Actual Unit Eigenvectors of A')
disp(V)
[ lambda,x ] = InversePM( A,x0,M )
end
