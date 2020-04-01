import numpy as np
def readsi(siname, sires=(24,24),sipts=1024, sisliceno=1):
    '''readsi - read siarray.1.1 file
    adapated from LoadSInD.m
    https://stackoverflow.com/questions/2146031/what-is-the-equivalent-of-fread-from-matlab-in-python
    https://stackoverflow.com/questions/44335749/read-a-float-binary-file-into-2d-arrays-in-python-and-matlab
    '''
    with open(siname, 'r') as fp1:
        # if not slice 1, calculate offset
        pxs = sires[0] * sires[1]
        offsetptr = pxs*sipts*4*2*(sisliceno-1)
        fp1.seek(offsetptr);
        # little-endian float32
        SI = np.fromfile(fp1, '<4f').reshape(pxs, 2*sipts).T
    return(SI)

def integrateSI(SI, sires, s, e):
     return np.sum(SI[s:e,:],0).reshape(sires)




