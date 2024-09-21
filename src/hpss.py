import numpy as np
def hpss(Y,M_dash=5,L_dash=5,iter=10):
    if (Y.ndim==1):
        return np.zeros((1,Y.shape[0]))
    elif (Y.shape[0] < 5) | (Y.shape[1] < 5):
        return np.zeros_like(Y)
    H = Y/np.sqrt(2)
    H_bar = np.zeros_like(H)
    P = Y/np.sqrt(2)
    P_bar = np.zeros_like(P)
    L = H.shape[0]
    M = H.shape[1]
    for i in range(iter):
        for m in range(M):
            H_bar[:,m] = 1/M_dash*(np.sum(H[:,max(0,m-M_dash):min(M,m+M_dash+1)],axis=1)-H[:,m])
        for l in range(L):
            P_bar[l,:] = 1/L_dash*(np.sum(P[max(0,l-L_dash):min(L,l+L_dash+1),:],axis=0)-P[l,:])
        H = H_bar/np.sqrt(H_bar*H_bar+P_bar*P_bar+1e-5)*Y
        P = P_bar/np.sqrt(H_bar*H_bar+P_bar*P_bar+1e-5)*Y
    return P
