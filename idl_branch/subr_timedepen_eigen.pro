pro subr_timedepen_eigen,natom, eigenpath, nte, historyfile, output,$
  inifraction=inifraction,accuracy=accuracy, safe_fact=safe_fact,log=log
  ; Name:
  ;    timedepen_eigen: subroutine-time-dependent ionization calculations using eigenvalue method.
  ; Purpose:
  ;    This routine returns ion fractions for both non-equilibirum and equilibrium cases by solving time-dependent ionization
  ;    equations.
  ;
  ; (1) Input:
  ; natom:
  ;   Atomic number for calculation
  ; eigenpath:
  ;   The folder contianing eigen matrix
  ; nte:
  ;   The number of Te grids depending on eigen matrix, e.g., 501 or 5001.
  ; historyfile:
  ;   Temperature and density history information, which includs temperature, electron density and time serical:
  ;   Te: one dimensional array, in unit of Kelvin;
  ;   ne: one dimensional array, in unit of cm^-3;
  ;   time: one dimensional array, in unit of sencond.
  ; output:
  ;   Output file name
  ;
  ; (2) Parameters
  ; inifraction:
  ;   Initial ion fractions. The default is initial ionization equilibrium state.
  ; safe_fact:
  ;   Limit the time-step, e.g., 0.67
  ; accuracy:
  ;   Monitor the change of ion fractions during one time-step, e.g., 0.9999
  ; log:
  ;   If ture, then save all ion fraction at each time nodes.
  ;
  ; (3) Result
  ;   All ion fractions of this atom are saved in 'output.sav' file. It saves a IDL structure:
  ;   conce_strc = {nei:dblarr(30,np),$
  ;                 eqi:dblarr(30,np),$
  ;                 tpr:dblarr(np),$
  ;                 rho:dblarr(np)}.
  ;   Here nei means non-equilibrium case, eqi is for ionization equilibrium, and np is for index of different history.
  ;   For example, the ion fraction of FeXXI in NEI case for the first record could be read from conce_strc.nei[20, 0].
  ;
  ;
  ; Examples:
  ;
  ;
  ; Update:
  ;   Chengcai Shen,
  ;   2014-07-10
  ;   Save the fraction information at each time frame, and add pth_icharg.
  ;   2014-09-11
  ;   Input parameters:
  ;     nte
  ;     dlogte
  ;     te_table
  ;   2014-09-12
  ;   switch: sw_time_step_log for recording the real time nodes and total steps.
  
  ;--------------------------------------------------------------------
  ; Parameters
  ;--------------------------------------------------------------------
  if(not keyword_set(safe_fact))then begin
    safe_fact = 0.67
  endif
  if(not keyword_set(accuracy))then begin
    accuracy = 0.9999
  endif
  if(not keyword_set(inifraction))then begin
    ini_eqi = 1
  endif else begin
    ini_eqi = 0
    ; check inifraction
    if(abs(1.0d0-total(inifraction)) ge 1.0d-4)then begin
      print,'Check initial ion fractions.'
      stop
    endif
  endelse
  
  
  ;--------------------------------------------------------------------
  ;  Te table
  ;--------------------------------------------------------------------
  dlogte = 5.0/float(nte-1)
  te_table = 10.0^(4.0 + dindgen(nte)*dlogte)
  
  ; -------------------------------------------------------------------
  ;   (1) common block
  ; -------------------------------------------------------------------
  ;common com_inival,inii_com,initime_com,np_com
  elemt = func_elemt_iv(natom)
  matric_evals = dblarr(natom+1)
  matric_evect = dblarr(natom+1,natom+1)
  matric_evect_invers = dblarr(natom+1,natom+1)
  
  ; -------------------------------------------------------------------
  ;   (2) read eigenvalues and eigenvectors
  ; -------------------------------------------------------------------
  openr,lun,eigenpath+string(elemt)+'eigen.dat',/get_lun,/f77_unformatted
  ;readu,lun,nte
  ;te_arr      = dblarr(nte)
  eqistate    = dblarr(natom+1,nte)
  eigenvalues = dblarr(natom+1,nte)
  eigenvector = dblarr(natom+1,natom+1,nte)
  eigenvector_invers = dblarr(natom+1,natom+1,nte)
  c_rate      = dblarr(natom+1,nte)
  r_rate      = dblarr(natom+1,nte)
  ;readu,lun,te_arr
  readu,lun,eqistate
  readu,lun,eigenvalues
  readu,lun,eigenvector
  readu,lun,eigenvector_invers
  readu,lun,c_rate
  readu,lun,r_rate
  free_lun,lun
  
  ; -------------------------------------------------------------------
  ;   (3) main cycle, for each trajectory
  ; -------------------------------------------------------------------
  np = 1l
  ntime = 1l
  openr,lun,historyfile,/get_lun
  readf,lun,np
  
  for ip = 0, np - 1 do begin
  
    ;     read density and temperature history
    readf,lun,ntime
    time_array = dblarr(ntime)
    te_array   = dblarr(ntime)
    rho_array  = dblarr(ntime)
    readf,lun,te_array
    readf,lun,rho_array
    readf,lun,time_array
    
    ; define variables
    conce_strc = {nei:dblarr(30,np),$
                  eqi:dblarr(30,np),$
                  te:dblarr(np),$
                  rho:dblarr(np)}
    ;
    
    dt_est_tol = dblarr(ntime)
    dt_est_ne = dblarr(ntime)
    dt_est_te = dblarr(ntime)
    
    ; save log file
    if(keyword_set(log))then begin
      openw,lun_log,string(ip,'(I06)')+'.log',/get_lun
    endif
    
    ;--------------------------------------------------------------------------
    ; (A) equilibrium state at initial time, and the final time
    ;--------------------------------------------------------------------------
    
    ; final time
    index_te = func_index_nearest(te_table, te_array[ntime-1])
    conce_eqi = eqistate[0:natom, index_te]
    
    ; initial time
    if(ini_eqi eq 1)then begin
      index_te = func_index_nearest(te_table, te_array[0])
      conce_ini = eqistate[0:natom, index_te]
    endif else begin
      conce_ini[0:natom] = inifraction[0:natom]
    endelse
    
    ;--------------------------------------------------------------------------
    ; (B) non-equilibrium
    ;--------------------------------------------------------------------------
    ; 1. estimate time-step using eigenvalues and density
    ;         check ne*t & dne
    for i = 1-1, ntime-1 do begin
      index_te = func_index_nearest(te_table, te_array[i])
      rho = rho_array[i]
      eval_max = max(abs(eigenvalues(*,index_te)))
      dt_est_elem = accuracy/(eval_max * rho)
      dt_est_ne[i] = safe_fact*dt_est_elem
    endfor
    
    
    ; 2. estimate time-step using tempreture
    for i = 0, ntime - 2 do begin
      tec = te_array[i]
      te1 = te_array[i+1]
      dte = abs(alog10(te1) - alog10(tec))
      
      if(dte lt dlogte)then begin
        i_plus = 0
        dte = 0.0
        while (dte lt dlogte) and ((i+i_plus) le ntime-2) do begin
          i_plus = i_plus + 1
          te1 = te_array[i + i_plus]
          dte = abs(alog10(te1) - alog10(tec))
        endwhile
        dt_est_te[i] = time_array[i+i_plus] - time_array[i]
      endif else begin
        dtime = time_array[i+1] - time_array[i]
        dt_est_te[i] = dtime/(dte/dlogte)
      endelse
    endfor
    dt_est_te[0] = dt_est_te[1]
    dt_est_te[ntime-1] = dt_est_te[ntime-2]
    
    
    ; 3. dt_est
    for i = 0, ntime-1 do begin
      dt_est_tol[i] = max([dt_est_ne[i], dt_est_te[i]])
    endfor
    ;
    ; (3) enter main cycle
    
    ; the initial time along this "streamline"
    i_loc = 1 - 1
    time = time_array[i_loc]
    
    while (i_loc le ntime - 1 -1) do begin
      te_prev = te_array[i_loc]
      
      ; how many points can be skipped?
      i_plus = 0
      dt_est = min(dt_est_tol(i_loc: (i_loc + i_plus)))
      dt = 0.0
      while ((dt lt dt_est) and (i_loc + i_plus le (ntime - 2))) do begin
        i_plus = i_plus + 1
        dt = time_array[i_loc + i_plus] - time
        dt_est = min(dt_est_tol(i_loc: (i_loc + i_plus)))
      endwhile
      ;
      ; ! the number of skipped points
      i_plus = i_plus - 1
      if (i_plus le 0) then begin
        i_plus = 1
      endif
      if ((i_loc + i_plus) ge ntime-1) then begin
        i_plus = 1
      endif
      dt = time_array(i_loc + i_plus) - time
      
      ;
      ;--------------------------------------------------------------------------
      ; case 1: i_plus > 1, skip several points
      ;--------------------------------------------------------------------------
      if (i_plus ge 2)then begin
        ; !             estimate avarage of density
        rhot = 0.0d0
        for i = i_loc, i_loc + i_plus - 1 do begin
          rhot_temp = (time_array(i + 1) - time_array(i)) * 0.5d0 * (rho_array(i + 1) + rho_array(i))
          rhot = rhot + rhot_temp
        endfor
        rho_now = rhot/dt
        te_now = te_prev
        index_te = func_index_nearest(te_table, te_now)
        
        conce_nei = func_solveionization_eigen(dt,conce_ini,rho_now,eigenvalues[0:natom,index_te],$
          eigenvector[0:natom,0:natom,index_te],eigenvector_invers[0:natom,0:natom,index_te])
        conce_ini = conce_nei
        
      endif
      
      ;--------------------------------------------------------------------------
      ; case 2: i_plus = 1, interpolation
      ;--------------------------------------------------------------------------
      if (i_plus eq 1)then begin
        ;!!                 a. simple interpolation
        te_sta = alog10(te_prev)
        te_end = alog10(te_array(i_loc + 1))
        dd_te = te_end - te_sta
        n_inter = fix(abs(dd_te)/dlogte) + 3
        
        dd_time = dt/float(n_inter - 1)
        
        for i = 1, n_inter - 1 do begin
          time_inner = float(i - 1) * dd_time
          
          te_now = interpol([te_array(i_loc), te_array(i_loc + 1)],[0.0,dt],time_inner)
          rho_now = interpol([rho_array(i_loc), rho_array(i_loc + 1)],[0.0,dt],$
            time_inner + 0.5d0 * dd_time)
            
          index_te = func_index_nearest(te_table, te_now)
          conce_nei = func_solveionization_eigen(dd_time,conce_ini,rho_now,eigenvalues[0:natom,index_te],$
            eigenvector[0:natom,0:natom,index_te],eigenvector_invers[0:natom,0:natom,index_te])
          conce_ini = conce_nei
          
        endfor
      endif
      
      ;!             time add 'one/several step' when cycle 50/60 has been
      ;!             performed. It does not depend on the number of 'n_inter'.
      i_loc = i_loc + i_plus
      time = time_array(i_loc)
      
      ;    save log file
      if(keyword_set(log))then begin
        printf,lun_log,time
        printf,lun_log,conce_nei[0:natom]
      endif
      
    endwhile
    
    ;------------------------------------------------------------------
    ;    save the results
    ;------------------------------------------------------------------
    conce_strc.nei[0:natom,ip] = conce_nei[0:natom]
    conce_strc.eqi[0:natom,ip] = conce_eqi[0:natom]
    conce_strc.te[ip] = te_array[ntime-1]
    conce_strc.rho[ip] = rho_array[ntime-1]
    
    ;    print to screen
    curloop=strcompress(string(fix(ip+1),f='(i5)'),/rem)
    print,format='(%"\33[1m running time-dependent ionization (%s  %s / %d ).\33[1a")',string(elemt),curloop,np
    
  endfor ;  end ip cycle
  
  ;--------------------------------------------------------------------
  ;  save charge status of all streamline at the finnal time
  ;--------------------------------------------------------------------
  save,filename = output,conce_strc
  
  ;
  print,'sub_timedepen completed.'
end
