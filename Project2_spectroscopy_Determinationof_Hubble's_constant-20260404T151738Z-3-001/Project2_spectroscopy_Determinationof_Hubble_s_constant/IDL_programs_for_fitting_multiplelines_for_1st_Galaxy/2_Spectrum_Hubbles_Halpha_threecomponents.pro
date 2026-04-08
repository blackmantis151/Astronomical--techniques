filename = 'test.fits'

data = mrdfits(filename, 1, header)

;help, data, /structure
;IDL> help, data, /structure
;** Structure <27a4d848>, 8 tags, length=32, data length=32, refs=1:
;   FLUX            FLOAT           13.5624
;   LOGLAM          FLOAT           3.58120
;   IVAR            FLOAT         0.0988757
;   AND_MASK        LONG                 0
;   OR_MASK         LONG                 0
;   WDISP           FLOAT           1.04876
;   SKY             FLOAT           11.6919
;   MODEL           FLOAT           8.11338

flux = data.FLUX
loglam = data.LOGLAM
wave = 10.0^loglam
model=data.MODEL

;C III = 1908.7 Å
;C IV = 1548.2 Å
;Mg II = 2802.7 Å
;Hα = 6562.7 Å
;Hβ = 4861.3 Å
;Hγ= 4340.4 Å
;Hδ = 4101.7 Å
;[O III] = 4958.9, 5006.8 Å
;[N II] = 6544.2, 6583.5 Å
;[S II] = 6715.8, 6732.8 Å
;O II = 3727.3 Å

lines=['Halpha','Hbeta','Hgamma','Hdelta','NII2',  'SII1', 'SII2',  'OII', 'CIII', 'Mg II']
lamb0=[6562.7, 4861.3, 4340.4, 4101.7,    6583.5, 6715.8, 6732.8, 3727.3, 1908.7, 2802.7]

H0standard=67.4;km/s/Mpc
distance=1597.44;Mpc
velc=3e+5;km/s

plamb=lamb0*(1.0+H0standard*distance/velc)

print,'Predited Lamb, Rest frame Lamb0'
forprint,lines,plamb,lamb0,plamb-lamb0,textout=2

;lines=['Halpha','Hbeta','Hgamma','Hdelta','NII2',  'SII1', 'SII2',  'OII', 'CIII', 'Mg II']
;lamb0=[6562.7, 4861.3, 4340.4, 4101.7,    6583.5, 6715.8, 6732.8, 3727.3, 1908.7, 2802.7]

;uniform background around Halpha line
condi_bkg=where(wave ge 8720 and wave le 8800)
mean_bkg_model=mean(model(condi_bkg))
print,'background flux',mean_bkg_model

;condi=where(wave ge plamb(0)-170 and wave le plamb(0)+170)
condi=where(wave ge  8936.15-180.0 and wave le  8936.15+190.0); 8936.15 is decided based on bestfit
Hareq_wave=wave(condi)
Hareq_flux=flux(condi)-mean_bkg_model
Hareq_model=model(condi)-mean_bkg_model

tek_color

set_plot,'ps'
!P.FONT=0
device,set_font='Helvetica-Bold'
device,/inches,xsize=8.5,ysize=5.5
device,filename='Spectrum_Halpha_multiGaussfit.ps',/color
;multiplot,[1,2];,ygap=0.00004

plotsym,0,0.005,/fill
plot, wave, flux-mean_bkg_model, psym=8, color=0, xtitle='Wavelength (Angstrom)', $
      ytitle='Flux', linestyle=0, thick=0.001, xrange=[plamb(0)-200,plamb(0)+250],/xstyle,yrange=[-1,20],/ystyle

oplot,wave, model-mean_bkg_model, sym=-3, color=2, thick=5

;Halpha
oplot,[plamb(0),plamb(0)],[2,17],linestyle=2,color=4,thick=3
xyouts,plamb(0),15,'Halpha',charsize=1.2,orientation=90

;NII2
oplot,[plamb(4),plamb(4)],[2,17],linestyle=2,color=4, thick=3
xyouts,plamb(4),12,'NII-2',charsize=1.2,orientation=90

;SII1
oplot,[plamb(5),plamb(5)],[2,17],linestyle=2,color=4, thick=3
xyouts,plamb(5),12,'SII-1',charsize=1.2,orientation=90

;SII2
oplot,[plamb(6),plamb(6)],[2,17],linestyle=2,color=4, thick=3
xyouts,plamb(6),10,'SII-2',charsize=1.2,orientation=90

;;----------------
;Hafit = gaussfit(Hareq_wave,Hareq_model,Hafitparms,NTERMS=3)
;print,'Gauss fit parameters Halpha:'
;print,Hafitparms
;oplot,Hareq_wave,Hafit,psym=-3,color=3,thick=3
;-----------------

;;;8,12;9 heights
;;;8910,8930,8960; central wavelengths
;;;5, 10, 7; dispersons
ppt_start=[3.0,8910.0,5.0,   8.0,8930.0,10.0,   4.0,8960.0,7.0]
;;;

;central component
;ppt_start=[10.0,8930.0,10.0]

;ppt = replicate({value:0.D, fixed:0, limited:[0,0], $
;                limits:[0.D,0]}, 3)

ppt = replicate({value:0.D, fixed:0, limited:[0,0], $
                limits:[0.D,0]}, 9)

ppt[0].fixed=0
ppt[0].limited(0) =1
ppt[0].limits(0)=0.1;ppt_start(0)-5.00d
ppt[0].limits(1)=6.0;ppt_start(0)+5.00d

ppt[1].fixed=0
ppt[1].limited(0) =1
ppt[1].limits(0)=8900.0;ppt_start(1)-0.50D
ppt[1].limits(1)=8920.0;ppt_start(1)+0.750D

ppt[2].fixed=0
ppt[2].limited(0) =1
ppt[2].limits(0)=2.0;ppt_start(2)-0.5d
ppt[2].limits(1)=7.0;ppt_start(2)+0.750d

;-------------
;ppt_start=[3.0,8910.0,5.0,   8.0,8930.0,10.0,   4.0,8960.0,7.0]

ppt[3].fixed=0
ppt[3].limited(0) =1
ppt[3].limits(0)=2.0;ppt_start(0)-5.00d
ppt[3].limits(1)=12.0;ppt_start(0)+5.00d
;;
ppt[4].fixed=0
ppt[4].limited(0) =1
ppt[4].limits(0)=8910.0;ppt_start(1)-0.50D
ppt[4].limits(1)=8950.0;ppt_start(1)+0.750D
;;
ppt[5].fixed=0
ppt[5].limited(0) =1
ppt[5].limits(0)=3.0;ppt_start(2)-0.5d
ppt[5].limits(1)=12.0;ppt_start(2)+0.750d


;ppt_start=[3.0,8910.0,5.0,   8.0,8930.0,10.0,   4.0,8960.0,7.0]
ppt[6].fixed=0
ppt[6].limited(0) =1
ppt[6].limits(0)=0.1;ppt_start(0)-5.00d
ppt[6].limits(1)=7.0;ppt_start(0)+5.00d

ppt[7].fixed=0
ppt[7].limited(0) =1
ppt[7].limits(0)=8940.0;ppt_start(1)-0.50D
ppt[7].limits(1)=8970.0;ppt_start(1)+0.750D

ppt[8].fixed=0
ppt[8].limited(0) =1
ppt[8].limits(0)=5.0;ppt_start(2)-0.5d
ppt[8].limits(1)=9.0;ppt_start(2)+0.750d

ppt[*].value=ppt_start;[40.0,2.2,1.0, 35.0,4.0,2.5]
pweights=fltarr(n_elements(hareq_wave))+1.0d
prerr=fltarr(n_elements(hareq_wave))
;;
pfitpara=mpfitfun('makethreemultigauss',hareq_wave,hareq_model,err=prerr,start_params=ppt_start,parinfo=ppt,perror=pfiterror,yfit=pyfit,weights=pweights,niter=100,nfree=nfree,covar=covar,status=status,bestnorm=bestnorm,/silent)
;;
print,'----P fit params---'
;;print,'---Peak, Centar, sigma---'
forprint,pfitpara,pfiterror,format='(2(F15.2,1X))'
;;;print,pfiterror,format='(9(F15.2,1X))'
;;
make_grid,min(hareq_wave)-2,max(hareq_wave)+2,2.0,pgridd
;;
pfitresult1=[pfitpara(0),pfitpara(1),pfitpara(2)]
makegauss,pfitresult1,pgridd,p_gx1,p_gy1
oplot,p_gx1,p_gy1,linestyle=1,thick=2
;;
pfitresult2=[pfitpara(3),pfitpara(4),pfitpara(5)]
makegauss,pfitresult2,pgridd,p_gx2,p_gy2
oplot,p_gx2,p_gy2,linestyle=1,thick=2

pfitresult3=[pfitpara(6),pfitpara(7),pfitpara(8)]
makegauss,pfitresult3,pgridd,p_gx3,p_gy3
oplot,p_gx3,p_gy3,linestyle=1,thick=2

oplot,p_gx3,p_gy1+p_gy2+p_gy3,psym=-3,linestyle=2,thick=3
device,/close

Ha_shifted=pfitpara(4)
Velocity=(3e+5)*(Ha_shifted-lamb0(0))/lamb0(0)
Ha_redshift=(Ha_shifted-lamb0(0))/lamb0(0)
print,'Ha-shifted A, Ha A, Hashifted-Ha A, Vel (km/s), Redshift '
print,Ha_shifted,lamb0(0),Ha_shifted-lamb0(0),Velocity, Ha_redshift

end
