
FUNCTION aia_make_eis_raster, files, windata, no_sat=no_sat, offset=offset, $
                              clean=clean, eis_map=eis_map, verbose=verbose


;+
; NAME:
;     AIA_MAKE_EIS_RASTER
;
; PURPOSE:
;     Takes a sequence of AIA images and generates an EIS-style raster
;     image that matches the input EIS windata structure.
;
; CATEGORY:
;     SDO; AIA; Hinode; EIS; coalignment.
;
; CALLING SEQUENCE:
;     Result = AIA_MAKE_EIS_RASTER( Files, Windata )
;
; INPUTS:
;     Files:  List of AIA filenames (better to use cutouts than full-disk
;             images). It is strongly recommended that the maps be larger
;             than the EIS field-of-view and extend over a longer time
;             series than the EIS raster duration.
;     Windata:  An EIS data structure in the format produced by the
;               routine eis_getwindata.
;
; OPTIONAL INPUTS:
;     Offset:  Two-element array that specifies a spatial offset to apply
;              to the EIS data. That is, an EIS coordinate becomes
;              [x,y]+offset.
;     Eis_Map: An IDL map formed from an EIS line. The output map will be
;              formed from this map, only the image will be replaced by
;              the AIA synthetic raster. If not specified, then the
;              routine creates a map using the central wavelength of the
;              windata structure.
;	
; KEYWORD PARAMETERS:
;     NO_SAT:  Removes saturated images from the AIA sequence.
;     CLEAN:   Cleans the AIA images of cosmic rays (useful for "respiked"
;              AIA data).
;     VERBOSE: If set, then information about how many AIA images
;              contribute to each exposure is printed to the IDL window.
;
; CALLS:
;     SDO2MAP, EIS_MAKE_IMAGE, EIS_FIND_FILE, EIS_AIA_OFFSETS
;
; OUTPUTS:
;     Returns an IDL map that is the same size as one formed from the EIS
;     dataset, and has exactly the same coordinates. The image, however,
;     is formed from the AIA data. Note that the EIS map is corrected for
;     the offset given by eis_aia_offset.pro. To adjust the output map back
;     to the AIA coordinates, you need to add the OFFSET array (if this was
;     input).
;
; MODIFICATION HISTORY:
;     Ver.1, 27-Jan-2023, Peter Young
;     Ver.2, 25-Jul-2025, Peter Young
;       Fixed bug when the top of the EIS slit extends beyond the top of the
;       AIA image
;     Ver.3, 01-Aug-2025, Peter Young
;       Major change such that the alignment is now done on the sub-pixel
;       level by making use of inter_map.
;-


IF n_params() LT 2 THEN BEGIN
  print,'Use:  IDL> map=aia_make_eis_raster( files, windata [, /no_sat, offset= '
  print,'                   /clean, eis_map=, /verbose )'
  return,-1
ENDIF 

IF n_elements(offset) NE 2 THEN offset=[0.,0.]

;
; Load AIA data into maps.
;
map=sdo2map(files,no_sat=no_sat,clean=clean,/quiet)

aia_tai=anytim2tai(map.time)

;
; Get start and end times of each EIS exposure.
;
start_tai=anytim2tai(windata.time_ccsds)
END_tai=start_tai+windata.exposure_time

nx=windata.nx
ny=windata.ny

IF windata.hdr.slit_ind EQ 0 THEN swid=1.0 ELSE swid=2.0

;
; Get EIS coordinates.
;
xy=eis_aia_offsets(windata.hdr.date_obs)
eis_x=windata.solar_x+xy[0]+offset[0]
eis_y=windata.solar_y+xy[1]+offset[1]

s=size(map[0].data,/dim)

;
; I convolve the AIA image with a 2D gaussian of FWHM 3"
; 
nk=13
kernel=fltarr(nk,nk)
k=2.35482
ident=make_array(nk,value=1.)
kx=(findgen(nk)-nk/2)#ident
ky=ident#(findgen(nk)-nk/2)
r=sqrt(kx^2+ky^2)
kernel=exp(-r^2/2./(3.0/0.6/k)^2)

n_maps=n_elements(map)
FOR i=0,n_maps-1 DO BEGIN
  image=map[i].data
  newimage=convol(image,kernel)
  map[i].data=newimage
ENDFOR 

;
; Create the EIS map with eis_make_image if it wasn't input.
;
IF n_tags(eis_map) EQ 0 THEN BEGIN
  file=eis_find_file(windata.hdr_date_obs,/lev)
  eis_make_image,file,mean(windata.wvl),eis_map,/map
ENDIF

;
; This map has the dimensions of the EIS raster in the x direction,
; and the dimensions of the AIA image in the y-direction. It will
; contain the synthetic AIA raster, which will later by interpolated
; onto the EIS y-coordinates.
;
aia_eis_map=make_map(fltarr(nx,s[1]), $
                     xc=eis_map.xc, $
                     yc=map[n_maps/2].yc, $
                     dx=eis_map.dx, $
                     dy=map[n_maps/2].dy, $
                     id='AIA synthetic raster', $
                     time=eis_map.time, $
                     xunits=map[n_maps/2].xunits, $
                     yunits=map[n_maps/2].yunits, $
                     dur=map[n_maps/2].dur, $
                     roll_angle=map[n_maps/2].roll_angle, $
                     roll_center=map[n_maps/2].roll_center, $
                     l0=map[n_maps/2].l0, $
                     b0=map[n_maps/2].b0, $
                     rsun=map[n_maps/2].rsun)


;
; The following loops over the EIS exposures (i) and populates the
; columns of aia_eis_map.
; Typically there will be more than one AIA image that contributes
; to the exposure, and also more than column within the image.
; These are all averaged.
;
FOR i=0,nx-1 DO BEGIN
  ;
  ; Find AIA images that match the EIS exposure duration. I extend
  ; the duration by 2s in each direction to help ensure that I get
  ; some images.
  ;
  k=where(aia_tai GE start_tai[i]-2. AND aia_tai LE END_tai[i]+2.,nk)
  ;
  IF nk NE 0 THEN BEGIN
    count=0
    FOR j=0,nk-1 DO BEGIN
      ;
      ; Extract data columns from the AIA images that fall within
      ; the EIS slit width. I extend the width by 0.5 arcsec in both
      ; directions to help with the averaging.
      ;
      get_map_coord,map[k[j]],xp,yp
      ix=where(xp[*,0] GE eis_x[i]-swid/2.-0.5 AND xp[*,0] LE eis_x[i]+swid/2.+0.5,nix)
      IF nix EQ 0 THEN CONTINUE
      getmin=min(abs(eis_x[i]-xp),imin)
     ;
     ; Average over the columns and add to the output map.
     ;
      aia_eis_map.data[i,*]=aia_eis_map.data[i,*]+total(map[k[j]].data[ix,*],1)
      count=count+nix
    ENDFOR
    ;
    ; Normalize by the number of images used.
    ;
    IF count NE 0 THEN aia_eis_map.data[i,*]=aia_eis_map.data[i,*]/float(count)
    ;
    IF keyword_set(verbose) THEN print,format='(" Exposure: ",i3," -- maps: ",i2," -- columns: ",i3)',i,nk,count
    count=0
  ENDIF 
ENDFOR 

;
; Create a copy of eis_map and apply the y-offset to it.
;
eis_map2=eis_map
eis_map2.yc=eis_map2.yc+offset[1]

;
; Interpolate the AIA map onto the coordinate grid of the corrected EIS map.
;
new_map=inter_map(aia_eis_map,eis_map2)

return,new_map


;
; Get the 1D array of y-coordinates and add the OFFSET value.
; Note that yp is the same for all of the AIA maps, so I just
; take the most recent value (from the earlier i-loop).
; iy is the AIA pixel that corresponds to the bottom pixel of
; the EIS map.
;
;; get_map_coord,eis_map,exp,eyp
;; eyp=eyp+offset[1]
;; getmin=min(abs(eyp[0]-reform(yp[0,*])),iy)
;; print,getmin,iy

;; eis_map2=eis_map
;; eis_map2.yc=eis_map2

;
; Copy the EIS map into the final output map, and then replace the
; data with the synthetic EIS image from AIA.
;
;; outmap=eis_map
;; s=size(img,/dim)
;; IF s[1]-iy GE ny THEN BEGIN
;;   outmap.data=img[*,iy:iy+ny-1]
;; ENDIF ELSE BEGIN
;;   outmap.data[*,0:s[1]-iy-1]=img[*,iy:*]
;; ENDELSE

;; outmap.id='AIA raster'

;return,outmap

END
