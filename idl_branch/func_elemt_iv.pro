function func_elemt_iv,i
  case i of 
   0: str='0' 
   1: str='H'
   2: str='He'
   3: str='Li'
   4: str='Be'
   5: str='B'
   6: str='C'
   7: str='N'
   8: str='O'
   9: str='F'
   10: str='Ne'
   11: str='Na'
   12: str='Mg'
   13: str='Al'
   14: str='Si'
   15: str='P'
   16: str='S'
   17: str='Cl'
   18: str='Ar'
   19: str='K'
   20: str='Ca'
   21: str='Sc'
   22: str='Ti'
   23: str='V'
   24: str='Cr'
   25: str='Mn'
   26: str='Fe'
   27: str='Co'
   28: str='Ni'
   29: str='Cu'
   30: str='Zn'
  endcase 
  return,str
end