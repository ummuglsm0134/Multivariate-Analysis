# Multivariate-Analysis
#Mean,Correlation.Covariance
options ls=75 ps=65 nodate nonumber;  
 
/* read mineral.dat file */ 
data bones; 
infile "/home/u59165872/MINERAL/mineral.dat"; 
input subject dradius radius dhumerus humerus dulna ulna; 
run; 
 
/* Using Proc IML compute the mean vector x-bar */ 
 
/* Create the matrix x */ 
proc iml; 
  use bones; 
  read all var _num_ into x [colname = varnames]; 
  n=nrow(x); p=ncol(x); 
 
/* Compute x-bar */ 
j=j(n,1,1); 
t_j=t(j); 
x_bar=(1/n)*t(j)*x; 
print x_bar; 
 
/* J and I matrices */ 
J1=j(n,n,1);  
I=I(n); 
 
/*Compute MLE, S, and R */ 
 
MLE_s=(1/n)*t(X)*(I-(1/n)*J1)*X; print MLE_s; 
s=(1/(n-1))*t(X)*(I-(1/n)*J1)*X; 
d=diag(s); 
sd=sqrt(d); 
inv_sd=inv(sd); 
CorrMatrix= inv_sd *s *inv_sd; print CorrMatrix; 
quit; 
run; 
 
/*Question 2 */ 
 
proc iml; 
s={9 4.2 12 2.4, 4.2 4 7 7.2, 12 7 25 2, 2.4 7.2 2 16}; 
n=nrow(s); p=ncol(s); 
d=diag(s); 
sd=sqrt(d); 
inv_sd=inv(sd); 
CorrMatrix= inv_sd *s *inv_sd; print CorrMatrix; 
quit; 
run;
