function makesinglegauss,GX,A,GY

;make_grid,-50,-2.0,0.05,GX

h1=A[0]
c1=A[1]
w1=A[2];/2.355; A[2] is fwhm; sigma=fwhm/2.355

Z1 = (GX-c1)/w1

GY1 = h1*EXP(-Z1^2.d/2.d)
;--------------------

GY=GY1

return, GY

end
