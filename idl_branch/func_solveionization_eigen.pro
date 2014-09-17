function func_solveionization_eigen,time,f0,electro_ne,evals,matric_evect,matric_evect_invers
  n = n_elements(f0)
  ft = dblarr(n)
  ;  eigenvalues matric, including ne * time
  matric_evals = dblarr(n,n)
  e = 2.71828182845904d0
  for i = 0,n-1 do begin
    matric_evals(i,i) = e^(evals(i)*time*electro_ne)
  endfor
  
  ;  matric operation
  matric_1 = MATRIX_MULTIPLY(matric_evals,matric_evect)
  matric_2 = MATRIX_MULTIPLY(matric_evect_invers,matric_1)

  ft = MATRIX_MULTIPLY(f0,matric_2)

  return,ft
end