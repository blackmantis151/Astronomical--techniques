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

;lines=['Halpha','Hbeta','Hgamma','Hdelta','NII2',  'SII1', 'SII2',  'OII', 'CIII', 'Mg II']
;lamb0=[6562.7, 4861.3, 4340.4, 4101.7,    6583.5, 6715.8, 6732.8, 3727.3, 1908.7, 2802.7]


;          0        1        2        3      4        5        6      7       8        9       10      11
lines=['Halpha','Hbeta','Hgamma','Hdelta','NII2',  'SII1', 'SII2',  'OII', 'OIII1', 'OIII2', 'NeIII', 'OI' ]
lamb0=[6562.7, 4861.3, 4340.4, 4101.7,    6583.5, 6715.8, 6732.8,   3727.3, 4959,     5007, 3869.0, 6300.0]

H0standard=67.4;km/s/Mpc
distance=1597.44;Mpc
velc=3e+5;km/s

plamb=lamb0*(1.0+H0standard*distance/velc)

print,'Predited Lamb, Rest frame Lamb0'
forprint,lines,plamb,lamb0,plamb-lamb0,textout=2

;lines=['Halpha','Hbeta','Hgamma','Hdelta','NII2',  'SII1', 'SII2',  'OII', 'CIII', 'Mg II']
;lamb0=[6562.7, 4861.3, 4340.4, 4101.7,    6583.5, 6715.8, 6732.8, 3727.3, 1908.7, 2802.7]

;since the entire spectrum within the mentioned range is noisy, we considered entire spectrum
condi_bkg=where(wave ge 6660.0 and wave le 6700.0)
mean_bkg_model=mean(model(condi_bkg))
print,'background flux',mean_bkg_model

;NeIII

;condi=where(wave ge plamb(0)-170 and wave le plamb(0)+170)
condi=where(wave ge 6810.0-100.0 and wave le  6810.0+100.0); 8936.15 is decided based on bestfit
Hareq_wave=wave(condi)
Hareq_flux=flux(condi)-mean_bkg_model
Hareq_model=model(condi)-mean_bkg_model

tek_color

set_plot,'ps'
!P.FONT=0
device,set_font='Helvetica-Bold'
device,/inches,xsize=8.5,ysize=5.5
device,filename='Spectrum_OIII2_zoom_Gaussfit.ps',/color
;multiplot,[1,2];,ygap=0.00004

plotsym,0,0.5,/fill
plot, wave, flux-mean_bkg_model, psym=-8, color=0, xtitle='Wavelength (Angstrom)', $
      ytitle='Flux', linestyle=0, thick=0.001, xrange=[6810-100,6810+100],/xstyle,yrange=[-3,10],/ystyle

oplot,wave, model-mean_bkg_model, sym=-3, color=2, thick=5


;plotsym,0,0.5,/fill
;plot, wave, flux, psym=8, color=0, xtitle='Wavelength (Angstrom)', $
;     ytitle='Flux', linestyle=0, thick=0.001, xrange=[5250-150,5250+150],/xstyle,yrange=[-1,20],/ystyle
;oplot,wave, model, sym=-3, color=2, thick=2

;OIII2
oplot,[plamb(9),plamb(9)],[2,10],linestyle=2,color=4,thick=3
xyouts,plamb(9),7,'OIII2',charsize=1.2,orientation=90

;----------------
;Hafit = gaussfit(Hareq_wave,Hareq_model,Hafitparms,NTERMS=3)
;print,'Gauss fit parameters Halpha:'
;print,Hafitparms
;oplot,Hareq_wave,Hafit,psym=-3,color=3,thick=3
;-----------------

;;;8,12;9 heights
;;;8910,8930,8960; central wavelengths
;;;5, 10, 7; dispersons
;;;ppt_start=[8.0,8910.0,5.0,   12.0,8930.0,10.0,   9.0,8960.0,7.0]
;;;

ppt_start=[3,6810,3.0]

ppt = replicate({value:0.D, fixed:0, limited:[0,0], $
                limits:[0.D,0]}, 3)

;;ppt = replicate({value:0.D, fixed:0, limited:[0,0], $
;;                limits:[0.D,0]}, 9)
;;
;;ppt[0].fixed=0
;;ppt[0].limited(0) =0
;;ppt[0].limits(0)=5.0;ppt_start(0)-5.00d
;;ppt[0].limits(1)=10.0;ppt_start(0)+5.00d
;;
;;ppt[1].fixed=0
;;ppt[1].limited(0) =0
;;ppt[1].limits(0)=8900.0;ppt_start(1)-0.50D
;;ppt[1].limits(1)=8920.0;ppt_start(1)+0.750D
;;
;;ppt[2].fixed=0
;;ppt[2].limited(0) =0
;;ppt[2].limits(0)=2.0;ppt_start(2)-0.5d
;;ppt[2].limits(1)=7.0;ppt_start(2)+0.750d
;;
;;;-------------
;;
;;;ppt_start=[8.0,8910.0,5.0,   12.0,8930.0,10.0,   9.0,8960.0,7.0]
;;
ppt[0].fixed=0
ppt[0].limited(0) =1
ppt[0].limits(0)=0.1;ppt_start(0)-5.00d
ppt[0].limits(1)=5.0;ppt_start(0)+5.00d
;;
ppt[1].fixed=0
ppt[1].limited(0) =1
ppt[1].limits(0)=6810.-10.;ppt_start(1)-0.50D
ppt[1].limits(1)=6810.+10;ppt_start(1)+0.750D
;;
ppt[2].fixed=0
ppt[2].limited(0) =1
ppt[2].limits(0)=0.5;ppt_start(2)-0.5d
ppt[2].limits(1)=5.0;ppt_start(2)+0.750d
;;
;;;---------------------
;;
;;;ppt_start=[8.0,8910.0,5.0,   12.0,8930.0,10.0,   9.0,8960.0,7.0]
;;
;;ppt[6].fixed=0
;;ppt[6].limited(0) =0
;;ppt[6].limits(0)=7.0;ppt_start(0)-5.00d
;;ppt[6].limits(1)=11.0;ppt_start(0)+5.00d
;;
;;ppt[7].fixed=0
;;ppt[7].limited(0) =0
;;ppt[7].limits(0)=8950.0;ppt_start(1)-0.50D
;;ppt[7].limits(1)=8970.0;ppt_start(1)+0.750D
;;
;;ppt[8].fixed=0
;;ppt[8].limited(0) =0
;;ppt[8].limits(0)=5.0;ppt_start(2)-0.5d
;;ppt[8].limits(1)=9.0;ppt_start(2)+0.750d
;;
;;;;-------------------
;;;
;;;Histogauss_P_PA_LDN1225_excluding_starsaroundNGC7538_and_NIRexcessstars.psppt[3].fixed=0
;;;ppt[3].limited(0) =1
;;;ppt[3].limits(0)=15.0;ppt_start(3)-5.00d
;;;ppt[3].limits(1)=35.0;ppt_start(3)+5.00d
;;;
ppt[*].value=ppt_start;[40.0,2.2,1.0, 35.0,4.0,2.5]
;;
pweights=fltarr(n_elements(hareq_wave))+1.0d
prerr=fltarr(n_elements(hareq_wave))
;;
pfitpara=mpfitfun('makesinglegauss',hareq_wave,hareq_model,err=prerr,start_params=ppt_start,parinfo=ppt,perror=pfiterror,yfit=pyfit,weights=pweights,niter=100,nfree=nfree,covar=covar,status=status,bestnorm=bestnorm,/silent)
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
oplot,p_gx1,p_gy1,linestyle=0,thick=4,color=4
;;
;;pfitresult2=[pfitpara(3),pfitpara(4),pfitpara(5)]
;;makegauss,pfitresult2,pgridd,p_gx2,p_gy2
;;oplot,p_gx2,p_gy2,linestyle=2
;;
;;pfitresult3=[pfitpara(6),pfitpara(7),pfitpara(8)]
;;makegauss,pfitresult3,pgridd,p_gx3,p_gy3
;;oplot,p_gx3,p_gy3,linestyle=2

device,/close

Ha_shifted=pfitpara(1)
Ha_redshift=(Ha_shifted-lamb0(9))/lamb0(9)
Velocity=(3e+5)*(Ha_redshift)

print,'For OIII2 line----'

print,Ha_shifted,lamb0(9),Ha_shifted-lamb0(9),Velocity, Ha_redshift

end
