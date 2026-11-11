"""Read selected native statistic columns/feature indices without giant arrays."""
import zipfile
import numpy as np

def _read_exact(f,n):
    chunks=[];remaining=n
    while remaining:
        b=f.read(min(remaining,8*1024*1024))
        if not b:raise EOFError('Truncated native NPY payload')
        chunks.append(b);remaining-=len(b)
    return b''.join(chunks)

def _discard(f,n):
    while n:
        b=f.read(min(n,8*1024*1024))
        if not b:raise EOFError('Truncated native NPY payload')
        n-=len(b)

def selected_native_columns(path,indices,names,expected_feature_sha,expected_feature_count):
    indices=np.asarray(indices,dtype=np.int64)
    with np.load(path,allow_pickle=False) as meta:
        cols=meta['columns'].tolist();actual=str(meta['feature_order_sha256'])
        if actual!=expected_feature_sha:raise ValueError('Native feature-order hash mismatch')
        if len(cols)!=len(set(cols)) or len(names)!=len(set(names)):
            raise ValueError('Duplicate declared or requested statistic columns')
        if not names:raise ValueError('No requested statistic columns')
    wanted=[cols.index(name) for name in names]
    with zipfile.ZipFile(path) as z, z.open('values.npy') as f:
        version=np.lib.format.read_magic(f);shape,fortran,dtype=np.lib.format._read_array_header(f,version)
        if len(shape)!=2 or dtype.kind!='f':raise ValueError('Unsupported native stats array')
        if shape[0]!=expected_feature_count:raise ValueError('Native feature count mismatch')
        if shape[1]!=len(cols):raise ValueError('Native statistic width mismatch')
        if len(indices) and (indices.min()<0 or indices.max()>=shape[0]):raise IndexError('Feature index outside native array')
        result=np.empty((len(indices),len(wanted)),dtype=dtype)
        if fortran:
            for col in range(max(wanted)+1):
                if col in wanted:
                    vec=np.frombuffer(_read_exact(f,shape[0]*dtype.itemsize),dtype=dtype)
                    result[:,wanted.index(col)]=vec[indices]
                else:_discard(f,shape[0]*dtype.itemsize)
        else:
            start=0
            while start<shape[0]:
                nr=min(65536,shape[0]-start);chunk=np.frombuffer(_read_exact(f,nr*shape[1]*dtype.itemsize),dtype=dtype).reshape(nr,shape[1]);mask=(indices>=start)&(indices<start+nr)
                result[mask]=chunk[indices[mask]-start][:,wanted];start+=nr
    return result
