pro makegauss,A,xgrid,GX,GY
;inputs A fit params [peak, center, sigma]
;xgrid=x values

h=A[0]
c=A[1]
w=A[2];/2.355; A[2] is fwhm; sigma=fwhm/2.355

;make_grid,-50,-2.0,0.05,GX

GXX=xgrid

Z = (GXX-c)/w
GY = h*EXP(-Z^2./2. )
GX=xgrid

end
