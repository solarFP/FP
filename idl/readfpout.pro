PRO readfpout,fpout, atmfile, file  =file
on_error,2
if (n_elements(file) eq 0) then file ='out.h5'
fpout = readhdf5(file)
makefpoutstr,fpout,atmfile
END
