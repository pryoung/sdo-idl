
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
; img is used to help create the synthetic EIS image. The x-size
; is set to the EIS map x-size. For the y-direction, we initially
; set this to the size of the AIA image in arcsec. For example,
; if the AIA image is 200 pixels, then we get the y-size as 120
; pixels (200*0.6).
;
img=fltarr(nx,round(s[1]*map[0].dx))

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
; I go through each EIS exposure and do the following procedure.
;  - get AIA maps that lie within start/end times of exposure.
;  - convolve each AIA map with the kernel
;  - rebin the AIA map to match EIS pixel sizes
;  - find AIA data column nearest to EIS slit position and add
;    to output image
;  - if multiple maps for an exposure, average the data columns
;
FOR i=0,nx-1 DO BEGIN
  k=where(aia_tai GE start_tai[i]-2. AND aia_tai LE END_tai[i]+2.,nk)
  IF nk NE 0 THEN BEGIN
    count=0
    FOR j=0,nk-1 DO BEGIN
     ;
     ; Rebin the AIA map to the EIS y-pixel size (1 arcsec)
     ;
      s=size(map[k[j]].data,/dim)
      r_ny=round(s[1]*map[k[j]].dy)
      rmap=rebin_map(map[k[j]],s[0],r_ny)

     ;
     ; Find data column(s) that matches EIS location
     ;
      get_map_coord,rmap,xp,yp
      ix=where(xp[*,0] GE eis_x[i]-swid/2.-0.5 AND xp[*,0] LE eis_x[i]+swid/2.+0.5,nix)
      IF nix EQ 0 THEN CONTINUE
      getmin=min(abs(eis_x[i]-xp),imin)

     ;
     ; Add column to output image.
     ;
      img[i,*]=img[i,*]+total(rmap.data[ix,*],1)
      count=count+nix

    ENDFOR
    IF count NE 0 THEN img[i,*]=img[i,*]/float(count)
    IF keyword_set(verbose) THEN print,format='(" Exposure: ",i3," -- maps: ",i2," -- columns: ",i3)',i,nk,count
    count=0
  ENDIF 
ENDFOR 

;
; Create the EIS map with eis_make_image if it wasn't input.
;
IF n_tags(eis_map) EQ 0 THEN BEGIN
  file=eis_find_file(windata.hdr_date_obs,/lev)
  eis_make_image,file,mean(windata.wvl),eis_map,/map
ENDIF

;
; Get the 1D array of y-coordinates and add the OFFSET value.
; Note that yp is the same for all of the AIA maps, so I just
; take the most recent value (from the earlier i-loop).
; iy is the AIA pixel that corresponds to the bottom pixel of
; the EIS map.
;
get_map_coord,eis_map,exp,eyp
eyp=eyp+offset[1]
getmin=min(abs(eyp[0]-reform(yp[0,*])),iy)

;
; Copy the EIS map into the final output map, and then replace the
; data with the synthetic EIS image from AIA.
;
outmap=eis_map
s=size(img,/dim)
IF s[1]-iy GE ny THEN BEGIN
  outmap.data=img[*,iy:iy+ny-1]
ENDIF ELSE BEGIN
ENDELSE


outmap.id='AIA raster'

return,outmap

END
