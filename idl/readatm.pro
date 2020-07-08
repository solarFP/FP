PRO readatm,file

if (n_elements(file) eq 0) then file = 'atm.dat'
nz = 0L
nions = 0L
nneutrals = 0L
openr,lu,file,/get_lun,/f77
readu,lu,nz,nions,nneutrals
zin = dblarr(nz)
tg = dblarr(nz)
bfield = dblarr(nz)
dni = dblarr(nions,nz)
dnn = dblarr(nneutrals,nz)
mions = dblarr(nions)
Zion = dblarr(nions)
Zn = dblarr(nneutrals)
Enion = dblarr(nneutrals)
readu,lu,zin,tg, bfield, dni, dnn, mions, Zion, Zn, Enion
free_lun,lu
; move everything up on level
(SCOPE_VARFETCH('zin', /ENTER, LEVEL=-1)) = zin
(SCOPE_VARFETCH('tg', /ENTER, LEVEL=-1)) = tg
(SCOPE_VARFETCH('bfield', /ENTER, LEVEL=-1)) = bfield
(SCOPE_VARFETCH('dni', /ENTER, LEVEL=-1)) = dni
(SCOPE_VARFETCH('dnn', /ENTER, LEVEL=-1)) = dnn
(SCOPE_VARFETCH('mions', /ENTER, LEVEL=-1)) = mions
(SCOPE_VARFETCH('Zion', /ENTER, LEVEL=-1)) = Zion
(SCOPE_VARFETCH('Zn', /ENTER, LEVEL=-1)) = Zn
(SCOPE_VARFETCH('Enion', /ENTER, LEVEL=-1)) = Enion
end
