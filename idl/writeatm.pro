PRO writeatm,zin,tg,bfield, dni, dnn, mions, Zion, Zn, Enion, file

if (n_elements(file) eq 0) then file = 'atm.dat'
nz = n_elements(zin)
nions = n_elements(mions)
nneutrals =n_elements(Zn)
openw,lu,file,/get_lun,/f77
writeu,lu,nz,nions,nneutrals
writeu,lu,zin,tg, bfield, dni, dnn, mions, Zion, Zn, Enion
free_lun,lu
end
