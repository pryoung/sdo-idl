
FUNCTION aia_rgb_table, wvl, reverse=reverse

;+
; NAME:
;      AIA_RGB_TABLE
;
; PURPOSE:
;      For IDL plot objects, it's necessary to specify a
;      'rgb_table'. This routine creates the rgb_table for the
;      specified wavelength.
;
; CATEGORY:
;      SDO/AIA, color tables, images.
;
; CALLING_SEQUENCE:
;      IDL> Result = AIA_RGB_TABLE( Wvl )
;
; INPUTS:
;      Wvl:    The wavelength (angstroms) for which the table is
;              required.
;
; KEYWORD PARAMETERS:
;      REVERSE: If set, then the color table is reversed.
;
; OUTPUTS:
;      The RGB table that can be used as input to, e.g., the image.pro
;      object of IDL.
;
; EXAMPLES:
;      IDL> rgb_table=aia_rgb_table(131)
;      IDL> p=plot_map_obj(map131, rgb_table=rgb_table)
;
; MODIFICATION HISTORY:
;      Ver.1, 21-Oct-2014, Peter Young
;      Ver.2, 6-Sep-2017, Peter Young
;        Added /reverse keyword.
;-


aia_lct,r,g,b,wavelnth=wvl,/load
rgb_table=[[r],[g],[b]]

IF keyword_set(reverse) THEN rgb_table=reverse(rgb_table,1)

return,rgb_table

END
