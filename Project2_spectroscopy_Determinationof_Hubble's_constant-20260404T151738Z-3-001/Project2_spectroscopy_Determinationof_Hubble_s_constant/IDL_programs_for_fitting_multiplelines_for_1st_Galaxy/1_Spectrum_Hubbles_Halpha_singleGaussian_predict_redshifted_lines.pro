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

tek_color

set_plot,'ps'
!P.FONT=0
device,set_font='Helvetica-Bold'
device,/inches,xsize=6,ysize=6
device,filename='Spectrum.ps',/color
;multiplot,[1,2];,ygap=0.00004

plotsym,0,0.001,/fill;,color=0
plot, wave, flux, xtitle='Wavelength (Angstrom)', $
      ytitle='Flux', linestyle=0,thick=0.001, xrange=[3600,9300],/xstyle,yrange=[0,20],/ystyle

oplot, wave, model, sym=-3, color=2, thick=2

device, /close

Ha0=6562.7
Ha=6800.0

;--------------to zoom a certain spectral feature

;set_plot, 'X

set_plot,'ps'
!P.FONT=0
device,set_font='Helvetica-Bold'
device,/inches,xsize=9.5,ysize=5.5
device,filename='Spectrum_zoom.ps',/color
;multiplot,[1,2];,ygap=0.00004

plotsym,0,0.001,/fill,color=1
plot, wave, flux, xtitle='Wavelength (Angstrom)', $
      ytitle='Flux', linestyle=0,thick=0.001, xrange=[3800,9500],/xstyle,yrange=[1,17],/ystyle


;plotsym,0,0.001,/fill;,color=0
;plot, wave, flux, xtitle='Wavelength (Angstrom)', $
;      ytitle='Flux', linestyle=0,thick=0.001, xrange=[6600,7000],/xstyle,yrange=[5,17],/ystyle

oplot, wave, model, sym=-3, color=2, thick=2

oplot,[Ha,Ha],[5,12],linestyle=2,color=4,thick=2

xyouts,Ha,09,'H-alpha',charsize=1.2,charthick=3,orientation=90

;Ha0=6562.7
;Ha=6745.0;6800.0

redshift=(Ha-Ha0)/Ha0

;print,redshift_Halpha

vel_Halpha=redshift*(3e+5);km/s

print,'---------------------------------------'
print,'Halpha--related'
print,'Halpha shift by .....:',Ha-Ha0
print,'Halpha redshift......:',redshift
print,'vel_Halpha......km/s.:',vel_Halpha
print,'--------------------------------------'


;NII0= 6544.2
NII0=6583.5

;Hb=4861.3
;Hb0=Hb+(1+(Ha0-Ha)/Ha)

NII=NII0*(1+redshift)

Hb0 = 4861.3
Hb=Hb0*(1+redshift)

Hg0=4340
Hg=Hg0*(1+redshift)

Hd0= 4101.7 
Hd=Hd0*(1+redshift)

print,'--------------------------------------'
;print,'NII related---'
;print,'RestNII wavelength....A:',NII0
;print,' ObsNII wavelength....A:',NII

print,'NII shifted by........A:',NII-NII0
print,'Halpha shifted by ...A:',Ha-Ha0
print,'Hbeta sifted by......A:',Hb-Hb0
print,'Hgamma shifted by ...A:',Hg-Hg0
print,'Hdelta shifted by ...A:',Hd-Hd0
print,'---------------------------------------'

;oplot,[NII,NII],[5,12],linestyle=2,color=4,thick=2
;xyouts,NII,9,'NII',charsize=1.2,charthick=3,orientation=90

oplot,[Hb,Hb],[5,12],linestyle=2,color=4,thick=2
xyouts,Hb,9,'Hbeta',charsize=1.2,charthick=3,orientation=90

oplot,[Hg,Hg],[5,12],linestyle=2,color=4,thick=2
xyouts,Hg,9,'Hgamma',charsize=1.2,charthick=3,orientation=90

oplot,[Hd,Hd],[5,12],linestyle=2,color=4,thick=2
xyouts,Hd,9,'Hdelta',charsize=1.2,charthick=3,orientation=90

device, /close


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

;[O III]	4959	4959 × 1.3606 ≈ 6746 Å
;[O III]	5007	5007 × 1.3606 ≈ 6811 Å

;Ne III] 3869 3869 × 1.3606 ≈ 5266 Å

;[O I] 6300 6300 × 1.3606 ≈ 8572 Å

lines=['Halpha','Hbeta','Hgamma','Hdelta','NII2',  'SII1', 'SII2',  'OII', 'OIII1', 'OIII2', 'NeIII', 'OI' ]
lamb0=[6562.7, 4861.3, 4340.4, 4101.7,    6583.5, 6715.8, 6732.8,   3727.3, 4959,     5007, 3869.0, 6300.0]

;------------------------

H0standard=67.4;km/s/Mpc
distance=1597.44;Mpc
velc=3e+5;km/s

plamb=lamb0*(1.0+H0standard*distance/velc)

print,'Predited Lamb, Rest frame Lamb0'
forprint,lines,plamb,lamb0,plamb-lamb0,textout=2

;LOADCT, 0
;TVLCT, 150, 150, 150, 100

tek_color

set_plot,'ps'
!P.FONT=0
device,set_font='Helvetica-Bold'
device,/inches,xsize=8.5,ysize=5.5
device,filename='Spectrum_predicted_redshifted_wavelengths.ps',/color
;multiplot,[1,2];,ygap=0.00004

plotsym,0,0.001,/fill
plot, wave, flux, psym=8, color=0, xtitle='Wavelength (Angstrom)', $
      ytitle='Flux', linestyle=0, thick=0.001, xrange=[3800,9300],/xstyle,yrange=[1,18],/ystyle

;plotsym,0,0.001,/fill;,color=0
;plot, wave, flux, xtitle='Wavelength (Angstrom)', $
;      ytitle='Flux', linestyle=0,thick=0.001, xrange=[6600,7000],/xstyle,yrange=[5,17],/ystyle
;LOADCT, 0
oplot, wave, model, sym=-3, color=2, thick=2

;Halpha
oplot,[plamb(0),plamb(0)],[2,17],linestyle=2,color=4,thick=3
xyouts,plamb(0),15,'Halpha',charsize=0.75,orientation=90

;Hbeta
oplot,[plamb(1),plamb(1)],[2,17],linestyle=2,color=4, thick=3
xyouts,plamb(1),15,'Hbeta',charsize=0.75,orientation=90

;Hgamma
oplot,[plamb(2),plamb(2)],[2,17],linestyle=2,color=4, thick=3
xyouts,plamb(2),15,'Hgamma',charsize=0.75,orientation=90

;Hdelta
oplot,[plamb(3),plamb(3)],[2,17],linestyle=2,color=4, thick=3
xyouts,plamb(3),15,'Hdelta',charsize=0.75,orientation=90

;NII2
oplot,[plamb(4),plamb(4)],[2,17],linestyle=2,color=4, thick=3
xyouts,plamb(4),12,'NII-2',charsize=0.75,orientation=90

;SII1
oplot,[plamb(5),plamb(5)],[2,17],linestyle=2,color=4, thick=3
xyouts,plamb(5),12,'SII-1',charsize=0.75,orientation=90

;SII2
oplot,[plamb(6),plamb(6)],[2,17],linestyle=2,color=4, thick=3
xyouts,plamb(6),10,'SII-2',charsize=0.75,orientation=90

;OII
oplot,[plamb(7),plamb(7)],[2,17],linestyle=2,color=4, thick=3
xyouts,plamb(7),10,'OII',charsize=0.75,orientation=90

;OIII1
oplot,[plamb(8),plamb(8)],[2,17],linestyle=2,color=4, thick=3
xyouts,plamb(8),10,'OIII1',charsize=0.75,orientation=90

;OIII2
oplot,[plamb(9),plamb(9)],[2,17],linestyle=2,color=4, thick=3
xyouts,plamb(9),10,'OIII2',charsize=0.75,orientation=90

;lines=['Halpha','Hbeta','Hgamma','Hdelta','NII2',  'SII1', 'SII2',  'OII', 'OIII1', 'OIII2', 'NeIII' ]
;lamb0=[6562.7, 4861.3, 4340.4, 4101.7,    6583.5, 6715.8, 6732.8,   3727.3, 4959,     5007, 3869.0]

;NeIII
oplot,[plamb(10),plamb(10)],[2,17],linestyle=2,color=4, thick=3
xyouts,plamb(10),10,'NeIII',charsize=0.75,orientation=90

;OI
oplot,[plamb(11),plamb(11)],[2,17],linestyle=2,color=4, thick=3
xyouts,plamb(11),10,'OI',charsize=0.75,orientation=90

device,/close

;lines=['Halpha','Hbeta','Hgamma','Hdelta','NII2',  'SII1', 'SII2',  'OII', 'CIII', 'Mg II']
;lamb0=[6562.7, 4861.3, 4340.4, 4101.7,    6583.5, 6715.8, 6732.8, 3727.3, 1908.7, 2802.7]

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
device,filename='Spectrum_Halpha_zoom_Gaussfit.ps',/color
;multiplot,[1,2];,ygap=0.00004

plotsym,0,0.5,/fill
plot, wave, flux-mean_bkg_model, psym=8, color=0, xtitle='Wavelength (Angstrom)', $
      ytitle='Flux', linestyle=0, thick=0.001, xrange=[plamb(0)-200,plamb(0)+250],/xstyle,yrange=[-1,20],/ystyle

oplot,wave, model-mean_bkg_model, sym=-3, color=2, thick=2

;Halpha
oplot,[plamb(0),plamb(0)],[2,17],linestyle=2,color=4,thick=3
xyouts,plamb(0),15,'Halpha',charsize=1.2,orientation=90

;;;NII2
;;oplot,[plamb(4),plamb(4)],[2,17],linestyle=2,color=4, thick=3
;;xyouts,plamb(4),12,'NII-2',charsize=1.2,orientation=90
;;
;;;SII1
;;oplot,[plamb(5),plamb(5)],[2,17],linestyle=2,color=4, thick=3
;;xyouts,plamb(5),12,'SII-1',charsize=1.2,orientation=90
;;
;;;SII2
;;oplot,[plamb(6),plamb(6)],[2,17],linestyle=2,color=4, thick=3
;;xyouts,plamb(6),10,'SII-2',charsize=1.2,orientation=90
;;
;;;OIII1
;;oplot,[plamb(8),plamb(8)],[2,17],linestyle=2,color=4, thick=3
;;xyouts,plamb(8),10,'OIII1',charsize=1.2,orientation=90
;;
;;;OIII2
;;oplot,[plamb(9),plamb(9)],[2,17],linestyle=2,color=4, thick=3
;;xyouts,plamb(9),10,'OIII2',charsize=1.2,orientation=90
;;
;----------------
Hafit = gaussfit(Hareq_wave,Hareq_model,Hafitparms,NTERMS=3)
print,'Gauss fit parameters Halpha:'
print,Hafitparms
oplot,Hareq_wave,Hafit,psym=-3,color=3,thick=3
;-----------------

;;;8,12;9 heights
;;;8910,8930,8960; central wavelengths
;;;5, 10, 7; dispersons
;;;ppt_start=[8.0,8910.0,5.0,   12.0,8930.0,10.0,   9.0,8960.0,7.0]
;;;

ppt_start=[10.0,8930.0,10.0]

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
ppt[0].limits(0)=5.0;ppt_start(0)-5.00d
ppt[0].limits(1)=14.0;ppt_start(0)+5.00d
;;
ppt[1].fixed=0
ppt[1].limited(0) =1
ppt[1].limits(0)=8920.0;ppt_start(1)-0.50D
ppt[1].limits(1)=8940.0;ppt_start(1)+0.750D
;;
ppt[2].fixed=0
ppt[2].limited(0) =1
ppt[2].limits(0)=7.0;ppt_start(2)-0.5d
ppt[2].limits(1)=12.0;ppt_start(2)+0.750d
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
oplot,p_gx1,p_gy1,linestyle=2,thick=3
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
Velocity=(3e+5)*(Ha_shifted-lamb0(0))/lamb0(0)
Ha_redshift=(Ha_shifted-lamb0(0))/lamb0(0)
print,'Ha-shifted A, Ha A, Hashifted-Ha A, Vel (km/s), Redshift '
print,Ha_shifted,lamb0(0),Ha_shifted-lamb0(0),Velocity, Ha_redshift

end
