FUNCTION parsehdf5str,r
  if (size(r,/type) ne 8) then return,-1
  for i = 0, n_tags(r) -1 do begin
    if (size(r.(i),/type) eq 8) then begin
      if (r.(i)._type eq 'GROUP') then begin
        val = parsehdf5str(r.(i))
        nme = r.(i)._name
        if (n_elements(thisstr) eq 0) then thisstr = create_struct(nme,val) else thisstr = create_struct(thisstr,nme,val)
      endif else if (r.(i)._type eq 'DATASET' or r.(i)._type eq 'ATTRIBUTE') then begin
        nme = r.(i)._name & val = r.(i)._data
        if (nme eq 'nE') then nme = 'nEn' ; renames nE since it is a keyword in IDL 
        if (n_elements(val) eq 1) then val = val[0]
        if (n_elements(thisstr) eq 0) then thisstr = create_struct(nme,val) else thisstr = create_struct(thisstr,nme,val)
      endif
    endif
  endfor
  return, thisstr
END

FUNCTION readhdf5,file
;forward_function parsehdf5str
on_error,2
r = h5_parse(file,/read)
rstr = parsehdf5str(r)

return,rstr
END
