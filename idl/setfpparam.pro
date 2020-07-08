PRO setfpparam,p,def,pres,dores, use_cmd = use_cmd

if (n_elements(use_cmd) eq 0) then use_cmd = 0
if (use_cmd and n_elements(p) gt 0) then p = p[0] $
else if (dores and n_elements(pres) gt 0) then p = pres[0] $
else if (n_elements(p) gt 0) then p = p[0] $
else p = def
case (size(def,/type)) of 
  2: p = FIX(P)
  3: p = LONG(p)
  4: p = FLOAT(p)
  5: p = DOUBLE(p)
  7: p = STRING(p)
endcase     
END
