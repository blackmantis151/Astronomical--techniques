function makethreemultigauss,GX,A,GY

;make_grid,-50,-2.0,0.05,GX

h1=A[0]
c1=A[1]
w1=A[2];/2.355; A[2] is fwhm; sigma=fwhm/2.355

Z1 = (GX-c1)/w1

GY1 = h1*EXP(-Z1^2.d/2.d)
;--------------------

h2=A[3]
c2=A[4]
w2=A[5];/2.355; A[2] is fwhm; sigma=fwhm/2.355

Z2 = (GX-c2)/w2
GY2 = h2*EXP(-Z2^2.d/2.d)
;----------------

h3=A[6]
c3=A[7]
w3=A[8];/2.355; A[2] is fwhm; sigma=fwhm/2.355

Z3 = (GX-c3)/w3
GY3 = h3*EXP(-Z3^2.d/2.d)


GY=GY1+GY2+GY3

return, GY

end
