import h5py, time, numpy as np
p="/Users/saljh8/Dropbox/Collaborations/Shunya/TGEN-spatial-IPF/TGEN_IPF_Xenium.h5ad"
f=h5py.File(p,"r")
t=time.time(); a=f["X"][500000:550000,:]; dt=time.time()-t
print(f"50k-row slab read: {dt:.2f}s -> {a.nbytes/1e6:.0f} MB, {a.nbytes/1e6/dt:.0f} MB/s")
print(f"projected full-X pass: {dt*1630319/50000:.0f}s")
f.close()
