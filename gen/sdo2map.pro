
function sdo2map, filename, clean=clean, nonorm=nonorm, no_sat=no_sat, $
                  no_rot=no_rot, trange=trange, tmatch=tmatch, $
                  no_shared_lib=no_shared_lib, quiet=quiet, $
                  no_pointing_check=no_pointing_check, $
                  nsatpix=nsatpix

;+
; NAME:
;     SDO2MAP()
;
; CATEGORY:
;     SDO; maps.
;
; PURPOSE: 
;     This routine reads a SDO file(s) into an IDL map. The key
;     features are:
;
;     1. The data array is divided by the exposure time (except for
;        HMI data).
;     2. The time tag in the output structure is set to T_OBS and
;        converted to UTC. (HMI data have times given in TAI format by
;        default.) 
;     3. If the cutout images have been re-spiked, then the /CLEAN
;        keyword removes the spikes but without removing the data.
;
; INPUTS:
;     Filename:  The name of a SDO data file. Can be an array.
;
; OPTIONAL INPUTS:
;     Tmatch:   An array of times for which images are required. This
;               is useful, for example, if you have a set of HMI
;               images and you want only those AIA images that are
;               closest in time to the HMI images. Note that the
;               returned map array will have the same dimensions as
;               TMATCH.
;
;     Trange:   A 2-element string array specifying a time range for
;               which images are needed.
;
;     Nsatpix:  If /no_sat is set, then this controls the threshold
;               at which an image is considered saturated. The file
;               header contains the keyword nsatpix that gives the
;               number of saturated pixels. The default value is 500.
;
; KEYW0RDS:
;     CLEAN:  If set, then the routine 'AIA CLEAN_CUTOUT_SEQUENCE' is
;             used to remove cosmic rays. It is intended that this be
;             used if the data have previously been re-spiked with
;             AIA_RESPIKE. 
;
;     NONORM: If set, then the images are not normalized with the
;             exposure time.
;
;     NO_SAT: If set, then exposures with more than 100 saturated
;             pixels are removed from the sequence.
;
;     NO_ROT: By default the routine corrects the AIA roll angle using
;             the ROT_MAP routine. If /no_rot is given then the
;             correction is not done.
;     QUIET:  If set, then information messages are not printed.
;
;     NO_POINTING_CHECK: By default, the routine checks if the files'
;             pointing information is up-to-date. If not, then a
;             message is printed and you should re-download the file to
;             get the best pointing data. Setting this keyword stops
;             the pointing being checked.
;
; OUTPUTS:
;     An IDL map containing the data from FILENAME. If FILENAME is an
;     array, then the returned map will be an array of map
;     structures. If an error is found, then -1 is returned.
;
; CALLS:
;     AIA_CLEAN_CUTOUT_SEQUENCE
;
; HISTORY:
;     Ver.1, 29-Jul-2013, Peter Young
;     Ver.2, 18-Dec-2013, Peter Young
;        Now filters out any images with zero exposure time; added
;        /NONORM keyword.
;     Ver.3, 23-Jul-2014, Peter Young
;        Converted data to float to avoid problems when dividing by
;        exposure time; also added /no_sat keyword to remove badly
;        saturated images from the sequence.
;     Ver.4, 15-Aug-2014, Peter Young
;        Now checks if the stated exposure time (exptime) is
;        consistent with the planned exposure time (aimgshce). If not
;        then aimgshce is used as the exposure time. (This is based on
;        advice by Paul Boerner.)
;     Ver.5, 3-Apr-2015, Peter Young
;        Added /use_shared_lib for read_sdo call to speed it up for
;        compressed images.
;     Ver.6, 11-Jun-2015, Peter Young
;        I've switched off use_shared_lib for HMI.
;     Ver.7, 1-Sep-2015, Peter Young
;        Now prints help message if no inputs given.
;     Ver.8, 17-Dec-2015, Peter Young
;        Added TMATCH= and TRANGE= keywords.
;     Ver.9, 29-Mar-2016, Peter Young
;        Introuced /NO_SHARED_LIB to switch off the /use_shared_lib
;        option in the call to read_sdo. Sometimes it causes the data
;        to be corrupted.
;     Ver.10, 1-Apr-2017, Peter Young
;        Now check if we have network connection before doing master
;        pointing check.
;     Ver.11, 6-Jun-2019, Peter Young
;        Changed call to aia_clean_cutout_sequence_2 to
;        aia_clean_cutout_sequence, after the former routine was
;        renamed.
;     Ver.12, 09-Apr-2021, Peter Young
;        Added /silent in call to read_sdo.
;     Ver.13, 01-Oct-2021, Peter Young
;        I've moved the /clean algorithm before the for loop
;        instead of after.
;     Ver.14, 21-Jul-2022, Peter Young
;        Added /quiet in call to aia_clean_cutout_sequence; added /quiet
;        keyword.
;     Ver.15, 04-Dec-2024, Peter Young
;        Added /no_pointing_check keyword.
;     Ver.16, 25-Jul-2025, Peter Young
;        Modified how saturated images are dealt with, and introduced
;        nsatpix= input.
;     Ver.17, 01-Aug-2025, Peter Young
;        Redefine n after images have been selected with trange.
;-


IF n_params() LT 1 THEN BEGIN
  print,'Use:  map=sdo2map( filename [, /clean, /nonorm, /no_sat, /no_rot, trange=, '
  print,'                   tmatch=, /no_shared_lib, /quiet, /no_pointing_check, '
  print,'                   nsatpix= ])'
  return,-1
ENDIF 


list=file_search(filename,count=n)
IF list[0] EQ '' THEN BEGIN
  message,/info,/cont,'FILENAME not found. Returning...'
  return,-1
ENDIF 

;
; An image is considered saturated if it has more than nsatpix saturated pixels.
; The image header contains the number of saturated pixels.
;
IF keyword_set(no_sat) THEN BEGIN 
  IF n_elements(nsatpix) EQ 0 THEN nsatpix=500
ENDIF ELSE BEGIN
  nsatpix=0
ENDELSE

;
; Read only the index initially in order to do some filtering.
;
read_sdo,list,index

;
; Filter out saturated images if /no_sat set.
;
IF nsatpix GT 0 THEN BEGIN
  k=where(index.nsatpix LT nsatpix,nk)
  IF nk EQ 0 THEN BEGIN
    message,/info,/cont,'All images are saturated! Try switching off /no_sat, or change the saturation threshold with nsatpix=. The current value is '+trim(round(nsatpix))+'. Returning...'
    return,-1
  ENDIF
  ;
  IF nk LT n THEN BEGIN
    IF NOT keyword_set(quiet) THEN message,/info,/cont,'Number of images removed due to saturation: '+ $
                                           trim(n-nk)+'.'
    list=list[k]
    index=index[k]
    n=nk
  ENDIF 
ENDIF 

;
; Do a check for the master pointing. The aia2wcsmin is not
; documented, but I believe it modifies the input index structure to
; insert the updated pointing information.
;
; At present I only check if the master file is different for the
; updated structure.
;
; Only do check on AIA.
;
; Since the check requires an internet connection, then I make
; sure we're connected.
;
inst=strmid(index[0].instrume,0,3)
net_chck=have_network()
IF inst EQ 'AIA' AND net_chck EQ 1 AND NOT keyword_set(no_pointing_check) THEN BEGIN 
  chckindex=aia2wcsmin(index[0])
  IF NOT keyword_set(quiet) THEN print,'  Current MPO: ',index[0].mpo_rec
  IF NOT keyword_set(quiet) THEN print,'  Latest MPO:  ',chckindex.mpo_rec
  IF chckindex.mpo_rec NE index[0].mpo_rec THEN BEGIN
    print,'***WARNING: an updated master pointing is available for this data-set***'
  ENDIF 
ENDIF


;
; Filter index based on TRANGE= input.
;
IF n_elements(trange) NE 0 THEN BEGIN
  tobs_tai=anytim2tai(index.t_obs)
  tr0_tai=anytim2tai(trange[0])
  tr1_tai=anytim2tai(trange[1])
  k=where(tobs_tai GE tr0_tai AND tobs_tai LE tr1_tai,nk)
  index=index[k]
  list=list[k]
  n=nk
ENDIF


;
; If TMATCH has been specified, then reduce index and list to those
; elements that are closest to the times in TMATCH. Note that index
; and list will have the same size as TMATCH and it's possible
; that the same image may end up appearing twice in the output (if the
; cadence of TMATCH is higher than that of the SDO sequence).
;
IF n_elements(tmatch) NE 0 THEN BEGIN
  tmatch_tai=anytim2tai(tmatch)
  tobs_tai=anytim2tai(index.t_obs)
  n=n_elements(tmatch)
  ind=-1
  FOR i=0,n-1 DO BEGIN
    getmin=min(abs(tobs_tai-tmatch_tai[i]),imin)
    ind=[ind,imin]
  ENDFOR
  ind=ind[1:*]
  index=index[ind]
  list=list[ind]
ENDIF 

;
; Check if have a mixed data-set.
;
inst=strmid(index.instrume,0,3)
k1=where(inst EQ 'AIA',n1)
k2=where(inst EQ 'HMI',n2)
;
IF n1 NE 0 AND n2 NE 0 THEN BEGIN
  print,'% SDO2MAP: The file selection contains both AIA and HMI images. Please give a list'
  print,'           of only AIA images or only HMI images.'
  return,-1
ENDIF 

IF n1 NE 0 THEN inst='aia' ELSE inst='hmi'

;
; Filter out any images with zero exposure time.
;
IF tag_exist(index,'exptime') THEN BEGIN 
  k=where(index.exptime GT 0.,nk)
  IF nk NE 0 THEN list=list[k]
  IF nk NE n THEN BEGIN
    IF NOT keyword_set(quiet) THEN print,'% SDO2MAP:  '+trim(n-nk)+' files were removed due to zero exposure time.'
  ENDIF 
ENDIF 

;
; Now read the image files using the updated file list.
;
; PRY, 3-Apr-2015: I've added /use_shared_lib below as this
; speeds up reading of compressed files considerably.
;
; 11-Jun-2015: use_shared_lib seems to be mess up the HMI data,
; so I've switched it off in this case.
;
; 17-Dec-2015: I've switched use_shared_lib back on for HMI
;
; 29-Mar-2016: OK, it looks like use_shared_lib doesn't work
; for older HMI files, so I've introduced no_shared_lib to
; switch it off. 
;
IF keyword_set(no_shared_lib) THEN use_shared_lib=0 ELSE use_shared_lib=1
read_sdo,list,index,data,use_shared_lib=use_shared_lib,/silent
index2map,index,float(data),map



;
; This is the planned exposure time. Need to convert from milliseconds
; to seconds.
;
IF inst EQ 'aia' THEN aimgshce=index.aimgshce/1000.


n=n_elements(list)

;
; 'bad_sat' indicates which of the maps is badly saturated (0-no,
; 1-yes) and removes them from the sequence if /no_sat is set. (Only
; for AIA.)
;
;bad_sat=bytarr(n)

;
; Scale the intensity array by the exposure time, but only for AIA
; data. 
;
; For HMI it's necessary to specify an exposure time, so I just use 1s.
;
IF inst EQ 'aia' THEN BEGIN 
 ;
 ; Here I clean the map images with aia_clean_cutout_sequence. The
 ; commented line shows how to switch to trace_unspike_time.
 ;
  IF keyword_set(clean) THEN BEGIN
;    outdata=trace_unspike_time(index,data)
    aia_clean_cutout_sequence,map.data,outdata,/fill,/quiet
    map.data=temporary(outdata)
  ENDIF 
 ;
  FOR i=0,n-1 DO BEGIN
   ;
    ;; IF keyword_set(no_sat) THEN BEGIN
    ;;   k_sat=where(map[i].data GE 16000,n_sat)
    ;;   IF n_sat GE 100 THEN bad_sat[i]=1b
    ;; ENDIF 
   ;
    IF keyword_set(nonorm) THEN BEGIN
      map[i].data=map[i].data
    ENDIF ELSE BEGIN 
     ;
     ; The following checks for an exposure time anomaly between the actual
     ; and planned exposure time. If there's a big difference then
     ; the planned exposure time is correct (!).
     ;
      exptime=map[i].dur
      chck=abs(exptime-aimgshce[i])/aimgshce[i]*100.
      IF chck GE 10 THEN BEGIN
        exptime=aimgshce[i]
        IF NOT keyword_set(quiet) THEN BEGIN
          print,'% SDO2MAP: exposure time is not consistent with planned exposure time for image '+trim(j)+'. Using planned exposure time.'
          print,aimgshce[i],map[i].dur
        ENDIF 
      ENDIF 
     ;
     ; Normalize the exposure 
     ;
      map[i].data=map[i].data/exptime
      map[i].dur=1.0
    ENDELSE 
   ;
   ; The following corrects the roll angle. I checked that the only tags
   ; that rot_map modifies are 'data' and 'roll_angle'.
   ;
   ; I confirmed that the angle specified in the crota2 keyword is the
   ; one that ends up in the roll_angle tag, and that this angle is
   ; given in degrees and gives the clockwise angle from solar-north.
   ;
   ; I confirmed that the following commands all give the same result:
   ;    rot_map(map[i],roll_angle=0.)
   ;    rot_map(map[i],-map[i].roll_angle)
   ;    rot_map(map[i],-map[i].roll_angle,rcenter=[map.xc,map.yc])
   ;
   ; In addition, the result of the map rotation is very close to the
   ; image generated by doing:
   ;    rot(img,-map.roll_angle,/interp)
   ; which is the routine used by aia_prep (via ssw_reg)
   ;
    IF NOT keyword_set(no_rot) AND map[i].roll_angle NE 0. THEN BEGIN
      rmap=rot_map(map[i],-map[i].roll_angle,rcenter=[map.xc,map.yc])
      map[i].data=rmap.data
      map[i].roll_angle=rmap.roll_angle
    ENDIF 
  ENDFOR
ENDIF ELSE BEGIN
 ;
 ; For HMI the only things to do are set the exposure time (.dur) and
 ; rotate the image.
 ;
  map.dur=1.0
  FOR i=0,n-1 DO BEGIN 
    IF NOT keyword_set(no_rot) AND map[i].roll_angle NE 0. THEN BEGIN
      rmap=rot_map(map[i],roll=0)
      map[i]=rmap
    ENDIF 
  ENDFOR 
ENDELSE 

;; IF keyword_set(no_sat) THEN BEGIN
;;   k=where(bad_sat EQ 0,nk)
;;   IF nk NE 0 THEN BEGIN 
;;     map=map[k]
;;     index=index[k]
;;     IF NOT keyword_set(quiet) THEN print,'% SDO2MAP:  removed '+trim(nk)+' images that have bad saturation.'
;;   ENDIF 
;; ENDIF 


map.time=index.t_obs

date=anytim2utc(index.t_obs,/date,/ccsds)
time=anytim2utc(index.t_obs,/ccsds,/time,/trunc)

map.time=trim(date)+' '+trim(time)

return,map

END
