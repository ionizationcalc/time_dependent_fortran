;  Name:
;    pro_wrtie_ionizrecomb_rate
;  purpose:
;    write ionization rate and recombination rate into unformatted file
;  use:
;    chianti database
;    ssw
;  History:
;    The last version calculated rates using fortran routines 
;    'cal_ion_recomb_rates' and database file 'eii_files' (from Yuen-Kun Ko).
;    That is based on CHIANTI 6.
;  Update:
;    2013-11-06
;  !NOTICE!: This routine saves all non-zero recombination rates into 'r'.
;    Their position (or 'index') does NOT fix the ionization state or ionzation equations. 
;    2014-05-21
;    Add elements "H".
;
pro pro_write_ionizrecomb_rate
  ;  
  unform_file = 'ionrecomb_rate.dat'
  nte = 501
  te_arr = 10.0d0^(4.0d0+dindgen(nte)*0.01d0)
  ;  list of elements
  nelemt = 28
  elemt_arr = ['h','he','li','be','b','c','n','o','f','ne','na','mg','al','si','p','s','ci','ar','k','ca','sc','ti','v','cr','mn','fe','co','ni']
  natom_arr = indgen(28)+1
  c = dblarr(30,30,nte)
  r = dblarr(30,30,nte)
  ;  call chianti subroutines
  for ite = 0,nte-1 do begin
    te = te_arr(ite)
    print,ite,'Te=',te
    for ielemt = 0,nelemt-1 do begin
      z = natom_arr(ielemt)
      ; for ioniztion rate:
      for icharg = 0,z-1 do begin
        ion = icharg + 1
        char_numb = strmid(string(ion),0,1,/reverse_offset)
        if(ion ge 10)then char_numb = strmid(string(ion),1,2,/reverse_offset)
        gname = "'"+elemt_arr(ielemt)+'_'+char_numb+"'"
        c(icharg,ielemt,ite) = ioniz_rate(gname,te,z=z,ion=ion)
      endfor
      ;  for recombination rate
      for icharg = 1,z do begin
        ion = icharg + 1
        char_numb = strmid(string(ion),0,1,/reverse_offset)
        if(ion ge 10)then char_numb = strmid(string(ion),1,2,/reverse_offset)
        gname = "'"+elemt_arr(ielemt)+'_'+char_numb+"'"
        ;  !notice!: here save all non-zero recombination rates into 'r'.
        ;  Their position do NOT fix the ionization state or ionzation equations. 
        r(icharg-1, ielemt,ite) = recomb_rate(gname,te,z=z,ion=ion)
        ;print,'recomb',gname,te,z,ion
      endfor
    endfor
  endfor
  
  ;  write into file
  openw,lun,/get_lun,unform_file,/f77_unformatted
  writeu,lun,c
  writeu,lun,r
  free_lun,lun
  
  print,'normal stop'
end
