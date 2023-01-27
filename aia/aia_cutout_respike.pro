

PRO aia_cutout_respike, filelist, istart=istart, spikelist=spikelist, silent=silent

;+
; NAME
;     AIA_CUTOUT_RESPIKE
;
; PROJECT
;     SDO/AIA
;
; PURPOSE
;     Respikes a sequence of AIA cutouts, saving the new images as
;     fits files in the subdirectory 'respike'.
;
; INPUTS
;     FILELIST  A string array containing a list of AIA filenames. It
;               is assumed that the files are all in the same
;               directory.
;
; OPTIONAL INPUTS
;     ISTART    If specified, then the processing begins with index
;               ISTART instead of 0. Sometimes the re-spiking fails
;               because the respike file is not found. My solution is
;               just to restart the routine, but starting with the
;               next file in the sequence.
;
;     SPIKELIST A list of spikes files that correspond with the cutout
;               list. The routine checks the times of the spikes files
;               and makes sure that the time matches the cutout file
;               to within 3 seconds. It doesn't check the
;               wavelength though! 
;
; KEYWORD PARAMETERS:
;     SILENT:  If set, then text messages to the IDL window will be
;              stopped. 
;
; OUTPUTS
;     Creates a new sub-directory called 'respike' and saves the
;     re-spiked images as FITS files in this directory. The respike
;     directory is located under the directory containing the original
;     files.
;
; MODIFICATION HISTORY
;     Ver.1, 30-Apr-2014, Peter Young
;     Ver.2, 16-Mar-2016, Peter Young
;        added /use_shared_lib in call to read_sdo.
;     Ver.3, 19-Apr-2016, Peter Young
;        added SPIKELIST= keyword.
;     Ver.4, 28-Mar-2021, Peter Young
;        when checking if 'spikes' is in the filenames, I've
;        added file_basename.
;     Ver.5, 14-Apr-2021, Peter Young
;        added keyword /silent to stop the text output to the IDL
;        window. 
;-

IF n_params() EQ 0 THEN BEGIN
  print,'Use:  IDL> aia_cutout_respike, filelist [, istart=, spikelist=, /silent]'
  print,''
  print,'  spikelist= list of local spikes files corresponding to filelist'
  print,'  istart=    index of filelist from which to start'
  print,'  The re-spiked files will go into the subdirectory respike.'
  return 
ENDIF

;
; By mistake I sometimes includes spikes files in FILELIST, so I here
; I check for this.
;
chck=strpos(file_basename(filelist),'spikes')
k=where(chck GE 0,nk)
IF nk GT 0 THEN BEGIN
  print,'% AIA_CUTOUT_RESPIKE: your input file list includes spikes files. Please remove them. Returning...'
  return
ENDIF 

;
; Get the times for the spike files.
;
IF n_elements(spikelist) NE 0 THEN BEGIN
  read_sdo,spikelist,index,/use_shared_lib,silent=silent
  t_spike_tai=anytim2tai(index.t_obs)
ENDIF 

;
; I assume that all the files belong in the same directory, and I get
; the directory from the first file in the list.
;
dir=file_dirname(filelist[0])
outdir=concat_dir(dir,'respike')
info=file_info(outdir)
IF info.directory EQ 0 THEN BEGIN
  file_mkdir,outdir
  IF NOT keyword_set(silent) THEN print,'% AIA_CUTOUT_RESPIKE: created the directory ',outdir
ENDIF


n=n_elements(filelist)

IF n_elements(istart) EQ 0 THEN istart=0

count=0

FOR i=istart,n-1 DO BEGIN 
  IF NOT keyword_set(silent) THEN print,'% AIA_CUTOUT_RESPIKE: Processing file '+trim(i)+' of '+trim(n-1)

  read_sdo,filelist[i],index,data,/use_shared_lib,silent=silent
  tobs_tai=anytim2tai(index.t_obs)

 ;
 ; If the spikelist exists then I check to make sure the spike file is
 ; within 3 seconds of the cutout file.
 ;
  IF n_elements(spikelist) NE 0 THEN BEGIN
    getmin=min(abs(tobs_tai-t_spike_tai),imin)
    IF getmin LE 3 THEN BEGIN
      read_sdo,spikelist[imin],spindex,spdata,/use_shared_lib,silent=silent
      aia_respike,index,data,outindex,outdata,ispikedd=spdata
    ENDIF ELSE BEGIN
      aia_respike,index,data,outindex,outdata   ;,respike_err=respike_err
    ENDELSE 
    ;; read_sdo,spikelist[i],spindex,spdata,/use_shared_lib
    ;; tsp_tai=anytim2tai(spindex.t_obs)
    ;; IF abs(tobs_tai-tsp_tai) GE 3.0 THEN BEGIN
    ;;   print,'% AIA_CUTOUT_RESPIKE: the times of the spikelist file and the cutout file do no match!'
    ;;   print,'                      The spikes file will be downloaded from JSOC...'
    ;;   junk=temporary(spdata)
    ;; ENDIF 
    ;; aia_respike,index,data,outindex,outdata,ispikedd=spdata
  ENDIF ELSE BEGIN 
    aia_respike,index,data,outindex,outdata   ;,respike_err=respike_err
  ENDELSE
  
  ;; IF respike_err EQ 1 THEN BEGIN
  ;;   print,'%AIA_CUTOUT_RESPIKE: failed to respike image.'
  ;;   count=count+1
  ;; ENDIF ELSE BEGIN 
  aia_write_fits,outindex,outdata,outdir=outdir
  s=size(outdata,/dim)
  IF s[1] EQ 3 THEN stop
;  ENDELSE 
ENDFOR 

; The spike files are downloaded to the current working directory, so
; in the following I delete the downloaded files.
;
IF n_elements(spikelist) EQ 0 THEN BEGIN 
   IF NOT keyword_set(silent) THEN print,'% AIA_CUTOUT_RESPIKE: deleting the downloaded spikes files...'
   list=file_search('*spikes.fits')
   file_delete,list
   IF NOT keyword_set(silent) THEN print,'                     ...completed.'
ENDIF 

IF count NE 0 THEN BEGIN
  print,'% AIA_CUTOUT_RESPIKE: there were '+trim(count)+' images that could not be re-spiked.'
ENDIF 

END
