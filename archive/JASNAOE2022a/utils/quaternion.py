import numpy as np


def convert_z_axis(vec, z_axis):
    inv_norm = 1/np.linalg.norm(z_axis)
    mx = z_axis[0]*inv_norm
    my = z_axis[1]*inv_norm
    mz = z_axis[2]*inv_norm
    #   
    theta = np.arccos(mz)
    cos_theta = np.cos(theta * 0.5)
    sin_theta = np.sin(theta * 0.5)
    norm = np.sqrt(mx * mx + my * my)
    # 
    if norm > 0:
        qr = cos_theta
        qi = sin_theta * -my / norm
        qj = sin_theta * mx / norm
    else:
        qr = 1
        qi = 0
        qj = 0
    qri = qr * qi
    qij = qi * qj 
    qrj = qr * qj
    qxx = qr * qr + qi * qi - qj * qj
    qyy = qr * qr - qi * qi + qj * qj
    qzz = qr * qr - qi * qi - qj * qj
    # 
    x = vec[:, 0]
    y = vec[:, 1]
    z = vec[:, 2]
    # 
    xx = x * qxx + 2 * (y * qij + z * qrj)
    yy = y * qyy - 2 * (z * qri - x * qij)
    zz = z * qzz - 2 * (x * qrj - y * qri)
    # 
    sample = np.concatenate([
        xx[:, np.newaxis],
        yy[:, np.newaxis],
        zz[:, np.newaxis]],
        axis=1)
    return sample
