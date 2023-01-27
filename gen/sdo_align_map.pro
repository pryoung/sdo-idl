
FUNCTION sdo_align_map, map, xrange=xrange, yrange=yrange, isub=isub, $
                        roll_angle=roll_angle, tref=tref

;+
; NAME:
;     SDO_ALIGN_MAP()
;
; PURPOSE:
;     A sequence of SDO cutouts (such as downloaded from the JSOC) are
;     extracted from the original images without interpolation
;     performed. This means that if you make a movie, you will see
;     occasional 1-pixel jumps in the move. This routine performs
;     bilinear interpolation to the image sequence so that the movie
;     will be smooth.
;
;     In addition, the XRANGE and YRANGE keywords allow sub-regions to
;     be extracted from the interpolated sequence. This ensures that
;     any movie made from this sequence will also be smooth.
;
;     Note: this routine will work on any sequence of maps, not just
;     SDO. 
;
; CATEGORY:
;     SDO; image analysis.
;
; CALLING SEQUENCE:
;	Result = SDO_ALIGN_MAP( Map )
;
; INPUTS:
;     MAP:   An array of maps. It is recommended that this is created
;            with sdo2map.pro.
;
; OPTIONAL INPUTS:
;     XRANGE: A 2 element vector that specifies a sub-range in the
;             X-direction to be extracted. This is an alternative to
;             using the SUB_MAP routine.
;
;     YRANGE: A 2 element vector that specifies a sub-range in the
;             X-direction to be extracted. This is an alternative to
;             using the SUB_MAP routine.
;
;     ISUB:   Used in conjunction with XRANGE and YRANGE. By default
;             XRANGE and YRANGE are assumed to correspond to the first
;             frame in the sequence. ISUB is used to specify the index
;             of the frame to which XRANGE and YRANGE should apply.
;
;     TREF:   Equivalent of ISUB, but a time is input to specify the
;             reference frame.
;
;     ROLL_ANGLE:  This specifies a roll angle to be applied to the
;                  images, which can be useful if you want to align a
;                  solar feature to the X or Y axis. Note that the X
;                  and Y coordinates assigned to the map are no longer
;                  valid.  The roll is applied after  the alignment
;                  has been performed, but before XRANGE and/or YRANGE
;                  are implemented.
;
; OUTPUTS:
;     Returns a new map array structure containing the interpolated
;     data. If a problem is found, then -1 is returned.
;
; EXAMPLE:
;     Create an aligned map from a set of cutouts.
;    
;     IDL> list=file_search('*.fits')
;     IDL> map=sdo2map(list)
;     IDL> map_align=sdo_align_map(map)
;
; MODIFICATION HISTORY:
;     Ver.1, 18-Dec-2013, Peter Young
;     Ver.2, 27-May-2014, Peter Young
;        fixed bug when map contains only 1 frame.
;     Ver.3, 6-Apr-2015, Peter Young
;        updated header (no change to code); added ROLL_ANGLE=
;        optional input.
;     Ver.4, 16-Dec-2015, Peter Young
;        modified to work on a single map.
;     Ver.5, 18-Jan-2016, Peter Young
;        fixed bug for implementation of tref= that affected first
;        frame in movie.
;-


IF n_params() LT 1 THEN BEGIN
  print,''
  print,'Use:  IDL> newmap=sdo_align_map(map [, xrange=, yrange=, isub=, roll_angle= ])'
  print,''
  print,'  xrange, yrange - specified in arcsec for the first frame'
  print,'  roll_angle - apply a roll angle to images'
  print,'  isub   - allows a different reference frame to be specified for xrange, yrange'
  print,''
  return,-1
ENDIF 

n=n_elements(map)

outmap=map

;
; Identify the reference map, which is frame 0 unless tref is
; specified. 
;
tmap_tai=anytim2tai(map.time)
IF n_elements(tref) NE 0 THEN BEGIN
  tref_tai=anytim2tai(tref)
  getmin=min(abs(tmap_tai-tref_tai),isub)
  refmap=map[isub]
  tref_tai=anytim2tai(refmap.time)
ENDIF ELSE BEGIN
  isub=0
  refmap=map[isub]
  tref_tai=anytim2tai(refmap.time)
ENDELSE 
x0=refmap.xc
y0=refmap.yc

s=size(refmap.data,/dim)
nx=s[0]
ny=s[1]

IF n_elements(roll_angle) NE 0 THEN BEGIN
  roll_map=rot_map(refmap,roll=roll_angle)
ENDIF 

;
; The following creates the new, interpolated map sequence (OUTMAP).
;
FOR i=0,n-1 DO BEGIN
 ;
 ; Get center of current map
 ;
  x1=map[i].xc
  y1=map[i].yc
 ;
 ;  The following rotates the center of the reference map to the time
 ;  of the current map ( [x0,y0] -> xy )
 ;
  t1_tai=anytim2tai(map[i].time)
  dt_days=(t1_tai-tref_tai)/86400.
  lonlat=xy2lonlat([x0,y0],refmap.time)
  dlon=diff_rot(dt_days,lonlat[1])
  lonlat[0]=lonlat[0]+dlon
  xy=lonlat2xy(lonlat,map[i].time)
 ;
 ; dx,dy are the differences between the actual center of the map, and
 ; the rotated center of the reference map
 ;
  dx=x1-xy[0]
  dy=y1-xy[1]
 ;
 ; Get differences in pixel numbers
 ;
  dxpix=dx/map[i].dx
  dypix=dy/map[i].dy
 ;
 ; Create new, interpolated image for map.
 ;
  img=map[i].data
  ix=findgen(nx)-dxpix
  iy=findgen(ny)-dypix
  newimg=bilinear(img,ix,iy)
 ;
 ; Update the OUTMAP structure
 ;
  outmap[i].data=newimg
  outmap[i].xc=xy[0]
  outmap[i].yc=xy[1]
 ;
  IF n_elements(roll_angle) NE 0 THEN BEGIN
    rmap=rot_map(outmap[i],roll=roll_angle)
    roll_map=[roll_map,rmap]
  ENDIF 
ENDFOR

IF n_elements(roll_angle) NE 0 THEN outmap=temporary(roll_map)


;------
; Now extract sub_map if XRANGE and YRANGE are set.
;
; Since the previous code has aligned all the maps relative to each
; other, then it's just a case of extracting the same pixel
; region from each map. This region is [i0:i1] and [j0:j1].
;
IF n_elements(xrange) NE 0 OR n_elements(yrange) NE 0 THEN BEGIN
  dx=outmap[isub].dx
  dy=outmap[isub].dy
 ;
  xc=outmap[isub].xc
  yc=outmap[isub].yc
 ;
  IF n_elements(xrange) NE 0 THEN BEGIN
    x0=xrange[0]
    x1=xrange[1]
   ;
    i0=(x0-xc)/dx + (float(nx)-1)/2.
    i0=round(i0)
    IF i0 LT 0 THEN i0=0
   ;
    i1=(x1-xc)/dx + (float(nx)-1)/2.
    i1=round(i1)
    IF i1 GT nx-1 THEN i1=nx-1
  ENDIF ELSE BEGIN
    i0=0 & i1=nx-1
  ENDELSE 
 ;
  IF n_elements(yrange) NE 0 THEN BEGIN
    y0=yrange[0]
    y1=yrange[1]
   ;
    j0=(y0-yc)/dy + (float(ny)-1)/2.
    j0=round(j0)
    IF j0 LT 0 THEN j0=0
   ;
    j1=(y1-yc)/dy + (float(ny)-1)/2.
    j1=round(j1)
    IF j1 GT ny-1 THEN j1=ny-1
  ENDIF ELSE BEGIN
    j0=0 & j1=ny-1
  ENDELSE
 ;
 ; In the following I convert OUTMAP to an index and then create the
 ; new map (OUTMAP2) from this index and the new DATA array
 ;
  data=outmap.data[i0:i1,j0:j1]
;  data=data[i0:i1,j0:j1,*]
  map2index,outmap,index
  index2map,index,data,outmap2
 ;
 ; I now have to create new xcen and ycen values for each image.
 ;
  FOR i=0,n-1 DO BEGIN
    x0=(i0-(float(nx)-1)/2.)*dx + outmap2[i].xc
    x1=(i1-(float(nx)-1)/2.)*dx + outmap2[i].xc
    outmap2[i].xc=(x1+x0)/2.
   ;
    y0=(j0-(float(ny)-1)/2.)*dy + outmap2[i].yc
    y1=(j1-(float(ny)-1)/2.)*dy + outmap2[i].yc
    outmap2[i].yc=(y1+y0)/2.
  ENDFOR
 ;
 ; The following three tags are not copied over by map2index so they
 ; need to be manually copied here.
 ;
  outmap2.id=outmap.id
  outmap2.b0=outmap.b0
  outmap2.rsun=outmap.rsun
 ;
 ; Now replace outmap
 ;
  outmap=temporary(outmap2)
    
ENDIF 

return,outmap

END
