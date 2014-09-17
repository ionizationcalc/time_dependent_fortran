function func_index_nearest, arr, a0
  n = n_elements(arr)
  
  if(arr[0] le arr[n-1])then begin
    res = where(arr ge a0)
  endif else begin
    res = where(arr le a0)
  endelse
  
  if(res[0] eq -1)then begin
    index_fine = n-1
    return,index_fine
  endif else begin
    index = res[0]
  endelse
  
  iplus = -1
  dd1 = abs(arr(index) - a0)
  dd2 = abs(arr(index + iplus) - a0)
  
  if(dd1 le dd2)then begin
    index_fine = index
  endif else begin
    index_fine = index + iplus
  endelse
  return,index_fine
end
