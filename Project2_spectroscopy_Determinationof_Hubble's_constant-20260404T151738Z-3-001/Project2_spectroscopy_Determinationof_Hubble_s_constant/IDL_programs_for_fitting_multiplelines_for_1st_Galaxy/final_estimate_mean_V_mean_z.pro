readcol,'Data_Velocity_and_redshift.dat',f='(A,D,D,D,D,D)',line,lobs,l,shift,v,z

print,mean(v),median(v),stddev(v)
;       108074.09       108050.90       73.282637
print,mean(z),median(z),stddev(z)
;      0.36024695      0.36016967   0.00024427201
