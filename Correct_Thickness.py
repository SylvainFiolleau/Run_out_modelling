# This python script converts depths (DoDs) into thickness values (along normals)

import numpy as np


def Calc_Norm(source_dtm, source_dod, cellsize):
    ###### Normals calculation ######
    # calculates normals of the dtm based on the neighboors (E-W and N-S cells: b-d-f-h, ignoring the D8 corners: a-c-g-i)
    # inputs: dtm grid and its header information
    #    ___ ___ ___
    #   | a | b | c |
    #   |___|___|___|
    #   | d | e | f |
    #   |___|___|___|
    #   | g | h | i |
    #   |___|___|___|
    b = np.array(source_dtm)
    d = np.array(source_dtm)
    f = np.array(source_dtm)
    h = np.array(source_dtm)
    h = np.array(source_dtm)
    b[1:, :] = b[0:-1, :]  # shift cells down
    d[:, 1:] = d[:, 0:-1]  # shift cells right
    f[:, 0:-1] = f[:, 1:]  # shift cells left
    h[0:-1, :] = h[1:, :]  # shift cells up

    # calculate eigenvectors (vi = null(A-lambda_i*Identity_matrix))
    dx = b - h
    dy = f - d
    sz = (b - source_dtm) ** 2 + (d - source_dtm) ** 2 + (source_dtm - f) ** 2 + (source_dtm - h) ** 2
    # clearvars b d f h

    # v2 = N':
    ny = -(2 * dx) / ((4 * dx ** 2 + 4 * dy ** 2 + sz ** 2 - 4 * sz + 4) ** 0.5 - sz + 2)
    nx = -(2 * dy) / ((4 * dx ** 2 + 4 * dy ** 2 + sz ** 2 - 4 * sz + 4) ** 0.5 - sz + 2)
    # clearvars dx dy sz
    nz = np.ones((np.size(nx, 0), np.size(nx, 1)))

    # scale based on cellsize
    h_ratio = 1 / cellsize
    nx *= h_ratio
    ny *= h_ratio

    # normalise
    norm_n = (nx ** 2 + ny ** 2 + nz ** 2) ** 0.5
    nx /= norm_n
    ny /= norm_n
    nz /= norm_n
    # clearvars norm_n

    # flip if not pointing up
    nx[nz < 0] *= -1
    ny[nz < 0] *= -1
    nz[nz < 0] *= -1

    orthogonal_thickness = source_dod / norm_n

    # end of normal function with outputs: nx, ny and nz
    return orthogonal_thickness
