; photometry.pro
; *****************************************************************************

      pro photometry, obsdate

; -----------------------------------------------------------------------------
; Initialization
; -----------------------------------------------------------------------------
      fwhm  = 10		                         ; for the centroid calculation
      apers = [1.,2.,3.,4.,5.,6.,7.,8.,9.,10] ; aperture
      sky1  = 20.	                            ; sky annulus r_in
      sky2  = 40	                            ; sky annulus r_out
      apidx = 3                               ; Choose which aperture to use
      
; -------------------------------------------------------------------------------
; Construct file list in the 'obsdate' library
; -------------------------------------------------------------------------------
      PathData = './'+obsdate+'/'
      cd, PathData
      FileList = [file_search('IMG*.fits')] 
      cd, '..'
      nfile    = n_elements(FileList)	

; -----------------------------------------------------------------------------
; Read databases
; -----------------------------------------------------------------------------
      readcol,'../Databases/USF4observedtargets.txt',ra_ttau,dec_ttau,fibre, $
         cfg_field_name1,cfg_field_name2,obj_name,format='D,D,I,A,A,A'

      readcol,'../Databases/USF4guides.txt',source_id_guide,ra_guide,dec_guide, $
         phot_g_mean_mag_1,phot_bp_mean_mag_1,phot_rp_mean_mag_1,_2MASS,Jmag, $
         e_Jmag,Hmag,e_Hmag,Kmag,e_Kmag,objID,gmag,e_gmag,rmag,e_rmag,imag, $
         e_imag,zmag,e_zmag,ymag,e_ymag, $
         format='A,D,D,F,F,F,A,F,F,F,F,F,F,A,F,F,F,F,F,F,F,F,F,F'

      restore,'../Databases/apass.sav'

; -------------------------------------------------------------------------------
; Loop over frames
; -------------------------------------------------------------------------------
      for i=0,nfile-1 do begin
   
         print 
         print,'Working on ',strmid(FileList[i],0,19)
         print
         res    = readfits(PathData+FileList[i],h,/silent)
         filter = strtrim(sxpar(h,'FILTER'),2)
         mjd    = sxpar(h,'MJD-OBS')
;
; Photometry of T Tauri stars
;
         adxy, h, ra_ttau, dec_ttau, xx, yy
         ii = where((xx ge apers[apidx] and xx le 1024-apers[apidx]) and $
                    (yy ge apers[apidx] and yy le 1024-apers[apidx]) and $
                    (res[xx,yy] eq res[xx,yy]),nii)
         if(nii gt 0) then begin
            id_ttau         = strarr(nii)
            x_ttau          = fltarr(nii)
            y_ttau          = fltarr(nii)
            instmag_ttau    = fltarr(nii)	
            instmagerr_ttau = fltarr(nii)	
            stdmag_ttau     = fltarr(nii)	
            stdmagerr_ttau  = fltarr(nii)	
            for j=0, nii-1 do begin       
               id_ttau[j]   = obj_name[ii[j]]
               cntrd, res, xx[ii[j]], yy[ii[j]], xcen, ycen, fwhm
               if(xcen EQ -1) then begin
                  xcen = xx[ii[j]]
                  ycen = yy[ii[j]]
               endif 
               x_ttau[j]    = xcen
               y_ttau[j]    = ycen
               aper, res, xcen, ycen, mag, magerr, sky, skyerr, 1, apers, [sky1,sky2], [-65536.,65536.], /nan, /exact
               instmag_ttau[j]    = mag[apidx]
               instmagerr_ttau[j] = magerr[apidx]
               stdmag_ttau[j]     = -1.
               stdmagerr_ttau[j]  = -1.
            endfor
         endif
         nttau = nii
;
; Photometry of guide stars
;
         adxy, h, ra_guide, dec_guide, xx, yy
         ii = where((xx ge apers[apidx] and xx le 1024-apers[apidx]) and $
                    (yy ge apers[apidx] and yy le 1024-apers[apidx]) and $
                    (res[xx,yy] eq res[xx,yy]),nii)
         if(nii gt 0) then begin
            id_guide         = strarr(nii)
            x_guide          = fltarr(nii)
            y_guide          = fltarr(nii)
            instmag_guide    = fltarr(nii)	
            instmagerr_guide = fltarr(nii)	
            stdmag_guide     = fltarr(nii)	
            stdmagerr_guide  = fltarr(nii)	
            for j=0, nii-1 do begin       
               id_guide[j] = source_id_guide[ii[j]]
               cntrd, res, xx[ii[j]], yy[ii[j]], xcen, ycen, fwhm
               if(xcen EQ -1) then begin
                  xcen = xx[ii[j]]
                  ycen = yy[ii[j]]
               endif 
               x_guide[j]    = xcen
               y_guide[j]    = ycen
               aper, res, xcen, ycen, mag, magerr, sky, skyerr, 1, apers, [sky1,sky2], [-65536.,65536.], /nan, /exact
               instmag_guide[j]    = mag[apidx]
               instmagerr_guide[j] = magerr[apidx]
               if(filter eq 'g') then begin
                  stdmag_guide[j]    = gmag[ii[j]]
                  stdmagerr_guide[j] = e_gmag[ii[j]]
               endif
               if(filter eq 'r') then begin
                  stdmag_guide[j]    = rmag[ii[j]]
                  stdmagerr_guide[j] = e_rmag[ii[j]]
               endif
               if(filter eq 'i') then begin
                  stdmag_guide[j]    = imag[ii[j]]
                  stdmagerr_guide[j] = e_imag[ii[j]]
               endif
               if(filter eq 'z') then begin
                  stdmag_guide[j]    = zmag[ii[j]]
                  stdmagerr_guide[j] = e_zmag[ii[j]]
               endif
            endfor
         endif
         nguide = nii
;
; Photometry of APASS9 stars
;
         adxy, h, apass.raj2000, apass.dej2000, xx, yy
         ii = where((xx ge apers[apidx] and xx le 1024-apers[apidx]) and $
                    (yy ge apers[apidx] and yy le 1024-apers[apidx]) and $
                    (res[xx,yy] eq res[xx,yy]),nii)
         if(nii gt 0) then begin
            id_apass         = strarr(nii)
            x_apass          = fltarr(nii)
            y_apass          = fltarr(nii)
            instmag_apass    = fltarr(nii)	
            instmagerr_apass = fltarr(nii)	
            stdmag_apass     = fltarr(nii)	
            stdmagerr_apass  = fltarr(nii)	
            for j=0, nii-1 do begin    
               id_apass[j] = strtrim(string(apass[ii[j]].recno),2) 
               cntrd, res, xx[ii[j]], yy[ii[j]], xcen, ycen, fwhm
               if(xcen EQ -1) then begin
                  xcen = xx[ii[j]]
                  ycen = yy[ii[j]]
               endif 
               x_apass[j]    = xcen
               y_apass[j]    = ycen
               aper, res, xcen, ycen, mag, magerr, sky, skyerr, 1, apers, [sky1,sky2], [-65536.,65536.], /nan, /exact
               instmag_apass[j]    = mag[apidx]
               instmagerr_apass[j] = magerr[apidx]
               if(filter eq 'g') then begin
                  stdmag_apass[j]    = apass[ii[j]].g_mag
                  stdmagerr_apass[j] = apass[ii[j]].e_g_mag
               endif
               if(filter eq 'r') then begin
                  stdmag_apass[j]    = apass[ii[j]].r_mag
                  stdmagerr_apass[j] = apass[ii[j]].e_r_mag
               endif
               if(filter eq 'i') then begin
                  stdmag_apass[j]    = apass[ii[j]].i_mag
                  stdmagerr_apass[j] = apass[ii[j]].e_i_mag
               endif
            endfor
         endif            
         napass = nii
;
; Discard APASS stars which coincide with a T Tauri star
;             
         if(nttau ne 0 and napass ne 0) then begin
            iflag = intarr(napass) + 1
            for j=0,napass-1 do begin
               dd = sqrt((x_ttau-x_apass[j])^2+(y_ttau-y_apass[j])^2)
               if(min(dd) lt fwhm) then iflag[j] = 0
            endfor
            idx = where(iflag eq 1,nidx)
            if(nidx eq 0) then begin
               napass = 0
            endif else begin
               id_apass         = id_apass[idx]    
               x_apass          = x_apass[idx]  
               y_apass          = y_apass[idx]  
               instmag_apass    = instmag_apass[idx]  
               instmagerr_apass = instmagerr_apass[idx] 
               stdmag_apass     = stdmag_apass[idx]  
               stdmagerr_apass  = stdmagerr_apass[idx] 
               napass           = nidx
            endelse
         endif
;
; Plot image
;
         set_plot,'ps'
         PSFile = PathData + strmid(FileList[i],0,19) + '.eps'
         device,xsize=16.0,ysize=16.,file=PSFile, $
            /portrait,/Schoolbook,/color,bits=8,xoffset=3.0,yoffset=6,/isolatin1
         loadct,0
         plotsym,0,thick=7

         !p.position=[0,0,1,1]
         plot, [10000],[10000],xs=1,ys=1,xr=[-0.5,1024-0.5],yr=[-0.5,1024-0.5],$
         xthick=5,ythick=5,xcharsize=0.001,ycharsize=0.001,font=3

         zsort = res[where(res eq res)]
         zsort =zsort(sort(zsort))
         avgz = avg(zsort[0.05*n_elements(zsort):0.95*n_elements(zsort)])
         s = stddev(zsort[0.05*n_elements(zsort):0.95*n_elements(zsort)])
         tv, 255 - bytscl((res-avgz)+3*s, min=0, max=20*s)
         loadct, 4
         plot, [10000],[10000],xs=1,ys=1,xr=[-0.5,1024-0.5],yr=[-0.5,1024-0.5],$
           /noerase,xthick=5,ythick=5,xcharsize=0.001,ycharsize=0.001,font=3
         xyouts,0.5,1.03,FileList[i]+'   filter='+filter,charsize=1.7,charthick=5,$
            align=0.5,/normal,font=3
         loadct,4

         if(nttau ne 0) then begin
            oplot,[x_ttau],[y_ttau],psym=8,thick=5,symsize=2,col=50
            xyouts,x_ttau+20,y_ttau,id_ttau,font=3,charsize=0.7,col=50
         endif
         if(nguide ne 0) then begin
            oplot,[x_guide],[y_guide],psym=6,thick=5,symsize=1.5,col=100
            xyouts,x_guide+15,y_guide,id_guide,font=3,charsize=0.7,col=100
         endif
         if(napass ne 0) then oplot,[x_apass],[y_apass],psym=4,thick=3,symsize=1.5,col=150
         ;         
         device, /close
         PDFFile = strmid(PSFile,0,strlen(PSFile)-4) + '.pdf'
         spawn, 'ps2pdf ' + PSFile + ' ' + PDFFile
         comstring =  'pdfcrop --margins "5 5 5 5" --clip ' + PDFFile + ' ' + PDFFile
         spawn, comstring
         spawn, 'rm ' + PSFile
;
; Print the results in an ascii file
;             
         if(nttau ne 0 or nguide ne 0 or napass ne 0) then begin
            openw,1, PathData + strmid(FileList[i],0,19) + '_photo.txt'
               printf,1,'====================================================================================='
               printf,1,'                               Photometry of stars
               printf,1,'====================================================================================='
               printf,1,'Filter = '+filter,format='(A45)'
               printf,1,'Flag          ID            Type       x       y     Stdmag  Stderr  Instmag  Insterr'
               
               if(nttau ne 0) then begin
                  for j=0,nttau-1 do printf,1,'1',id_ttau[j],'  TTau  ',x_ttau[j],y_ttau[j], $
                     stdmag_ttau[j],stdmagerr_ttau[j],instmag_ttau[j],instmagerr_ttau[j], $
                     format='(A2,A23,A9,F8.2,F8.2,F9.3,F8.3,F9.3,F8.3)'
               endif
               if(nguide ne 0) then begin
                  for j=0,nguide-1 do printf,1,'1',id_guide[j],'  Guide ',x_guide[j],y_guide[j], $
                     stdmag_guide[j],stdmagerr_guide[j],instmag_guide[j],instmagerr_guide[j], $
                     format='(A2,A23,A9,F8.2,F8.2,F9.3,F8.3,F9.3,F8.3)'
               endif
               if(napass ne 0) then begin
                  for j=0,napass-1 do printf,1,'1',id_apass[j],'  APASS ',x_apass[j],y_apass[j], $
                     stdmag_apass[j],stdmagerr_apass[j],instmag_apass[j],instmagerr_apass[j], $
                     format='(A2,A23,A9,F8.2,F8.2,F9.3,F8.3,F9.3,F8.3)'
               endif
            close,1
         endif
;
; End of the loop
;          
      endfor

; -----------------------------------------------------------------------------
; The end
; -----------------------------------------------------------------------------
      stop
      end

; *****************************************************************************

