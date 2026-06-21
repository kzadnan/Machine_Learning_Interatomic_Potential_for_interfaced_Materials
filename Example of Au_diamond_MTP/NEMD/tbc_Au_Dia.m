ev=1.06e-5+5.23e-6;
step=0.0005;
Q=0.5*(ev)*1/step*1e12*1.6E-19;

b=7.62337 ;
c=15.24673 ;
Area=b*c;
Delta_T=324.71-270.57;


TBC_Gan_aln=Q*1e20/(Area*Delta_T)