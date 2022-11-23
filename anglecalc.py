# Script that finds the waveplate angles to implement the SU(2) unitary U using the reciprocal polarization gadget 
# Written by Teodor Strömberg

import numpy as np

''' Paulis '''
X = np.array([[0, 1],[1, 0]])
Y = np.array([[0, -1j],[1j, 0]])
Z = np.array([[1, 0],[0, -1]])
I = np.array([[1,0],[0,1]])



''' Angle finding function '''
''' Finds the waveplate angles to implement the SU(2) unitary U using the reciprocal polarization gadget '''
def identity_gadget_angles(U):

    # Ensure SU(2)
    U = U/np.sqrt(np.linalg.det(U).astype(complex))
    Ux = blochrot([np.pi],'X').dot(U) # Exponentiated X to stay in SU(2)
    b,a = np.linalg.eig(Ux)
    idx = np.argmax(np.angle(b)) # Pick the positive eigenphase
    theta = -np.angle(b[idx])*2 # rotation angle
    v1 = a[:,idx] # initial direction of rotation axis
    dmat_1 = state2dmat(v1)

    # Rotate the axis of the unitary to X
    tz = np.real(np.trace(Z.dot(dmat_1)))
    tx = np.real(np.trace(X.dot(dmat_1)))
    th1 = np.arctan2(tz,tx)
    v2 = blochrot([th1],'Y').dot(v1)


    dmat_2 = state2dmat(v2)
    ty = np.real(np.trace(Y.dot(dmat_2)))
    tx = np.real(np.trace(X.dot(dmat_2)))
    th2 = -np.arctan2(ty,tx)
    v3 = blochrot([th2],'Z').dot(v2)

    # Now find the angles for the outer waveplates
    # They should implement the unitary:
    Uqh = blochrot([th1, th2, np.pi/2],'YZX')

    H = np.array([1,0])
    L = np.array([1,1j])/np.sqrt(2)
    Hp = np.linalg.inv(Uqh).dot(H)
    Lp = np.linalg.inv(Uqh).dot(L)

    dmat_Hp = state2dmat(Hp)
    dmat_Lp = state2dmat(Lp)
    th_q,th_h = map2HL(dmat_Hp,dmat_Lp)


    th_q = 2*th_h - th_q; # exchanging the order of the HWP and QWP from HQ to QH

    # reduce QHH to QH (last H angle is 22.5)
    th_q_fw = th_q + 90;
    th_h_fw = th_q - th_h + 22.5 - 90;

    # same for the backwards dir
    # reduce HHQ to HQ
    th_q_bw = -th_q + 90;
    th_h_bw = 22.5 + th_h - th_q - 90;

    theta = theta - np.pi # offset beause we started with Ux
    theta_hwp = theta / 4 * 180 / np.pi + 90 # convert to degrees and add offset

    return np.mod(th_q_fw,180), np.mod(th_h_fw,180), np.mod(theta_hwp,180), np.mod(th_h_bw,180), np.mod(th_q_bw,180)

''' Constructs the identity gadget unitary from the waveplate angles '''
def identity_gadet(th_q_fw, th_h_fw, theta_hwp, th_h_bw, th_q_bw):

    # Positive and negative Faraday rotators:
    Fp = blochrot([np.pi/2],'Y')
    Fm = blochrot([-np.pi/2],'Y')

    # Outer waveplates on the left:
    Uqh_l = qwp(th_q_fw).dot(hwp(th_h_fw))

    # Outer waveplates on the right:
    Uqh_r = hwp(th_h_bw).dot(qwp(th_q_bw))

    # Asymmetric X-gadget:
    Uz = qwp(90).dot(hwp(theta_hwp)).dot(qwp(90))

    # Middle gadget (Faradays sandwiching the X):
    Um = Fm.dot(Uz).dot(Fp)

    # Everything together:
    U = Uqh_l.dot(Um).dot(Uqh_r)

    # Explicit contstruction
    # U = qwp(th_q_fw).dot(hwp(th_h_fw)).dot(Fm).dot(qwp(90)).dot(hwp(theta_hwp)).dot(qwp(90)).dot(Fp).dot(hwp(th_h_bw)).dot(qwp(th_q_bw))

    return U


''' Helper functions '''

''' quarter-wave plate at an angle theta (in degrees) from the vertical axis '''
def qwp(theta):
    theta = theta * np.pi / 180
    return blochrot([theta*2, np.pi/2, -theta*2],'YZY')

''' Half-wave plate at an angle theta (in degrees) from the vertical axis '''
def hwp(theta):
    theta = theta * np.pi / 180
    return blochrot([theta*2, np.pi, -theta*2],'YZY')

''' Returns the density matrix for a state vector '''
def state2dmat(state):
    return np.outer(state.flatten(),np.conj(state.flatten()))

''' Returns angles such that hwp(hwp_th)*qwp(qw_th) * |in> = |H>,
    where dmat = |in><in|, and also ensures that hwp(hwp_th)*qwp(qw_th) |L'> = |L> '''
def map2HL(dmat,dmat_lp):
    qwp_th = 1/2 * np.arctan2(np.real(np.trace(X.dot(dmat_lp))),np.real(np.trace(Z.dot(dmat_lp)))) * 180 / np.pi + 45

    U1 = qwp(qwp_th)
    U2 = np.conj(U1.T)

    dmat_p = U1.dot(dmat.dot(U2))
    tx = np.real(np.trace(X.dot(dmat_p)))
    tz = np.real(np.trace(Z.dot(dmat_p)))

    hwp_th = np.arctan2(tx,tz) * 180 / np.pi / 4

    qwp_th = np.mod(qwp_th,180)
    hwp_th = np.mod(hwp_th,90)

    return qwp_th, hwp_th

''' Generate Haar-random unitary '''
def rand_unitary(dim):
    X = (np.random.rand(dim,dim) + 1j*np.random.rand(dim,dim))/np.sqrt(2)
    Q,R = np.linalg.qr(X)
    R = np.diag(np.sign(np.diag(R)))
    U = Q.dot(R)
    return U

''' Example usage: blochrot([th1,th2,th3],'XYZ')
    returns: e^(-1i*th1*X) * e^(-1i*th2*Y) * e^(-1i*th3*Z) '''
def blochrot(angle,axes):

    l = len(axes)
    la = len(angle)
    if l != la:
        raise Exception('Number of angles doesn\'t match the number of rotation axes')

    U = I

    for r in range(l):
        axis = axes[r]
        if axis == 'X':
            P = X
        elif axis == 'Y':
            P = Y
        elif axis == 'Z':
            P = Z
        else:
            return -1
        U = U.dot(np.cos(angle[r]/2)*I - 1j*np.sin(angle[r]/2)*P)

    return U

if __name__ == '__main__':
    np.set_printoptions(precision=3)

    # generate random unitary to test the script
    U = rand_unitary(2)
    U = U / np.sqrt(np.linalg.det(U).astype(complex)) # ensure SU(2)


    # find the waveplate angles
    t1,t2,t3,t4,t5 = identity_gadget_angles(U)

    # construct the unitary in the two propagation directions
    U_fw = identity_gadet(t1,t2,t3,t4,t5)
    U_bw = identity_gadet(-t5,-t4,-t3,-t2,-t1)

    print('Waveplate angles: \n',np.round([t1,t2,t3,t4,t5],3),'\n')
    print('Constructed forward unitary: \n',np.round(identity_gadet(t1,t2,t3,t4,t5),3),'\n')
    print('Constructed backwards unitary: \n',np.round(identity_gadet(-t5,-t4,-t3,-t2,-t1),3),'\n')
    print('Original unitary: \n',U/np.sqrt(np.linalg.det(U).astype(complex)),'\n')
    print('Norm of U_fw - U_gadget: \n',np.linalg.norm(U_fw-U),'\n')
    print('Norm of U_bw - U_gadget: \n',np.linalg.norm(U_bw-U),'\n')
