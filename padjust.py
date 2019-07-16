import numpy as np

def padjust(p,method):

    '''
    adjust p values
    modelled on p.adjust from R

    Translated by: Serafeim Loukas, serafeim.loukas@epfl.ch, Jul 2019.
    '''

    p = np.array(p)
    n = float(len(p))
    p = p.reshape(-1,)

    if method == 'BH':
        tmp = p*n ; p[tmp < 1.] = tmp[tmp < 1.].copy(); p[tmp >= 1.] = 1.

    elif method == 'fdr':
            i = np.array(range(int(n),0,-1))
            o = p.argsort()[::-1]
            ro = o.argsort()
            x = np.array([float(n) / i]).T * p[o].reshape(-1,1)
        
            # cumulative minimum
            b = np.zeros(x.shape)
            b[0,0] = x[0,0].copy()
            x_flat = np.reshape(x, [1,-1], order='F').ravel()
            b_flat = np.reshape(b, [1,-1], order='F').ravel()
            for k in range(1,len(x)):
                b_flat[k] = np.minimum(x_flat[k],b_flat[k-1])
            p = b_flat[ro].copy() # reorder   
            
    else:
            raise Exception('Method not implemented')
    
    p[p < 1.] = p[p < 1.].copy(); p[p >= 1.] = 1.
    
    return p