PRO writesattr,val,nme,id
  aspace_id = h5s_create_scalar()
  type_id = h5t_idl_create(val)
  attr_id = h5a_create(id, nme, type_id, aspace_id)
  h5a_write,attr_id, val
  h5a_close,attr_id
  h5t_close,type_id
  h5s_close,aspace_id
END
PRO writearray,arr,nme,id
  dims = size(arr,/dim)
  dspace_id = h5s_create_simple(dims)
  type_id = h5t_idl_create(arr)
  dataset_id = h5d_create(id, nme, type_id, dspace_id)
  h5d_write,dataset_id, arr
  h5d_close,dataset_id
  h5t_close,type_id
  h5s_close,dspace_id
END

PRO writefpout,fpout, outfile =outfile

if (n_elements(outfile) eq 0) then outfile ='out.h5'
file_id = h5f_create(outfile)
group_id = h5g_create(file_id,'inputparams')
writesattr,fpout.inputparams.nEn,'nE',group_id
writesattr,fpout.inputparams.nmu,'nmu',group_id
writesattr,fpout.inputparams.nz,'nz',group_id
writesattr,fpout.inputparams.maxiter,'maxiter',group_id
writesattr,fpout.inputparams.patype,'patype',group_id
writesattr,fpout.inputparams.inc_relativity,'inc_relativity',group_id
writesattr,fpout.inputparams.inc_CC,'inc_CC',group_id
writesattr,fpout.inputparams.inc_synchro,'inc_synchro',group_id
writesattr,fpout.inputparams.inc_magmirror,'inc_magmirror',group_id
writesattr,fpout.inputparams.inc_RC,'inc_RC',group_id
writesattr,fpout.inputparams.oneD,'oneD',group_id
writesattr,fpout.inputparams.reflecttop,'reflecttop',group_id
writesattr,fpout.inputparams.reflectbottom,'reflectbottom',group_id
writesattr,fpout.inputparams.writeoutput,'writeoutput',group_id
writesattr,fpout.inputparams.Emin,'Emin',group_id
writesattr,fpout.inputparams.Emax,'Emax',group_id
writesattr,fpout.inputparams.tolres,'tolres',group_id
writesattr,fpout.inputparams.toldiff,'toldiff',group_id
writesattr,fpout.inputparams.implicit_theta,'implicit_theta',group_id
writesattr,fpout.inputparams.mbeam,'mbeam',group_id
writesattr,fpout.inputparams.zbeam,'zbeam',group_id
writesattr,fpout.inputparams.Ecut,'Ecut',group_id
writesattr,fpout.inputparams.dlt,'dlt',group_id
writesattr,fpout.inputparams.eflux,'eflux',group_id
writesattr,fpout.inputparams.pasigma,'pasigma',group_id
writesattr,fpout.inputparams.resist_fact,'resist_fact',group_id
h5g_close,group_id

group_id = h5g_create(file_id, 'atm')
writesattr,fpout.atm.nneutral,'nNeutral',group_id
writesattr,fpout.atm.nion,'nIon',group_id
writearray,fpout.atm.zin,'zin',group_id
writearray,fpout.atm.tg,'tg',group_id
writearray,fpout.atm.bfield,'bfield',group_id
writearray,fpout.atm.dni,'dni',group_id
writearray,fpout.atm.dnn,'dnn',group_id
writearray,fpout.atm.mion,'mion',group_id
writearray,fpout.atm.zion,'Zion',group_id
writearray,fpout.atm.zn,'Zn',group_id
writearray,fpout.atm.Enion,'Enion',group_id
h5g_close,group_id

group_id = h5g_open(file_id,'/')
writearray,fpout.E,'E',group_id
writearray,fpout.mu,'mu',group_id
writearray,fpout.z,'z',group_id
writearray,fpout.esvol,'esvol',group_id

writearray,fpout.heatrate,'heatrate',group_id
writearray,fpout.momrate,'momrate',group_id
writearray,fpout.f,'f',group_id
h5g_close,group_id
h5f_close,file_id
END
