


PRO aia_clean_cutout_sequence, data, outdata, fill=fill, quiet=quiet

;+
; NAME:
;     AIA_CLEAN_CUTOUT_SEQUENCE
;
; PURPOSE:
;     This routine is intended to remove cosmic rays from "re-spiked"
;     AIA data, without removing real data.
;
; CATEGORY:
;     SDO; AIA; cosmic rays.
;
; CALLING SEQUENCE:
;     AIA_CLEAN_CUTOUT_SEQUENCE, Data, OutData
;
; INPUTS:
;     Data:   A 3D array of AIA images.
;
; KEYWORD PARAMETERS:
;     FILL:  If set, then flagged pixels are filled with a median of
;            surrounding pixels (5x5 pixel region). Otherwise they
;            given the value of -100.
;     QUIET: If set, then information messages are not printed.
;
; OUTPUTS:
;     A 3D array of same size as DATA containing the cleaned images.
;
; RESTRICTIONS:
;     The routine is intended to be used on AIA full-disk or cutout
;     images that have been "re-spiked". It probably won't do
;     anything for regular data.
;
;     If the images have varying exposure times (e.g., during flares),
;     then the routine should be applied to the normalized data.
;
;     It is recommended that you use this routine through sdo2map by
;     using the /clean keyword. 
;
; PROGRAMMING NOTES:
;     My criteria for flagging pixels are:
;
;     1. The intensity in a pixel must be greater than at the same
;        pixel in both the preceding and subsequent exposures.
;     2. The intensity in a pixel must be greater than the average
;        intensity of the preceding and subsequent exposures by at
;        least a factor of FACTOR1 (see code).
;     3. The intensity in a pixel must be greater than the median of
;        the surrounding 5x5 pixel block by CUTOFF DN (see code).
;
;     For the first and last frame FACTOR1 is replaced by FACTOR2,
;     which is more stringent.
;
; MODIFICATION HISTORY:
;     Ver.1, 06-Jun-2019, Peter Young
;     Ver.2, 21-Jul-2022, Peter Young
;       Added /quiet; added parameter check; modified print statement.
;-

IF n_params() LT 1 THEN BEGIN
  print,'Use:  IDL> aia_clean_cutout_sequence, data, outdata [, /fill, /quiet]'
  return
ENDIF


cutoff=5.0
factor1=1.7
factor2=2.5

missing=-100

siz=size(data)

n=siz[3]

outdata=data

;
; This deals with first and last exposures. A stricter criterion is
; applied (factor2 is used instead of factor1).
;
FOR i=0,1 DO BEGIN 
  IF i EQ 0 THEN BEGIN
    jj=0 & kk=1
  ENDIF ELSE BEGIN
    jj=n-1 & kk=n-2
  ENDELSE 
  img1=float(data[*,*,jj])
  img2=float(data[*,*,kk])
  medimg5=fmedian(img1,5,5)
  avg_t=img2
  k=where( img1 GT img2 $
           AND abs((img1-avg_t)/avg_t) GT factor2 $
           AND img1-medimg5 GT cutoff, nk)

  IF nk NE 0 THEN BEGIN 
    outimg=img1
    IF keyword_set(fill) THEN BEGIN
      outimg[k]=medimg5[k]
    ENDIF ELSE BEGIN
      outimg[k]=missing
    ENDELSE 
    IF NOT keyword_set(quiet) THEN print,format='("Image ",i4," -- flagged pixels: ",i10)',i,nk
    outdata[*,*,jj]=outimg
  ENDIF
ENDFOR 


;
; All other exposures
;
FOR i=1,n-2 DO BEGIN
  img1=float(data[*,*,i-1])
  img2=float(data[*,*,i])
  img3=float(data[*,*,i+1])
 ;
  medimg5=fmedian(img2,5,5)
 ;
  avg_t=float(img1+img3)/2.
 ;
  k=where( img2 GT img1 AND img2 GT img3  $
          AND abs((img2-avg_t)/avg_t) GT factor1 $
          AND img2-medimg5 GT cutoff, nk)
  IF nk NE 0 THEN BEGIN 
    outimg=img2
    IF keyword_set(fill) THEN BEGIN
      outimg[k]=medimg5[k]
    ENDIF ELSE BEGIN
      outimg[k]=missing
    ENDELSE 
    IF NOT keyword_set(quiet) THEN print,format='("Image ",i4," -- flagged pixels: ",i10)',i,nk
    outdata[*,*,i]=outimg
  ENDIF 
ENDFOR

END
