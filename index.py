import numpy as np
import math
np.set_printoptions(precision=3)


def injecao_P(k):
    somaP = 0
    # soma_adj = 0
    for i in range(nbar):
        somaP += Vs[k]*Vs[i]*(G[k, i]*math.cos(theta[k] -
                                               theta[i]) + B[k, i]*math.sin(theta[k]-theta[i]))
        # if (adj_matrix[k, i] == 1 and k!=i):
        #    soma_adj += P_calculado(Vs[k], Vs[i], theta[k], theta[i], g_matrix[k,i], b_matrix[k,i], tap_matrix[k,i], tap_matrix[i,k])
    return somaP


def injecao_Q(k):
    somaQ = 0
    # soma_adj = 0
    for i in range(nbar):
        if (k != i):
            aux = Vs[k]*Vs[i]*(G[k, i]*math.sin(theta[k] -
                                                theta[i]) - B[k, i]*math.cos(theta[k]-theta[i]))
            # print(aux)
            somaQ += aux
        # if (adj_matrix[k,i] == 1 and k!=i):
        #    soma_adj +=  Q_calculado(Vs[k], Vs[i], theta[k], theta[i], g_matrix[k,i], b_matrix[k,i], bsh_matrix[k,i], tap_matrix[k,i], tap_matrix[i,k])
    somaQ += -pow(Vs[k], 2)*B[k, k]
    # soma_adj += Vs[k]*shunt_bar[k]/100
    return somaQ  # , soma_adj



def P_calculado(v1, v2, theta1, theta2, g, b):
    return v1*v1*g - (v1 * v2 * (g * math.cos(theta1 - theta2) + b * math.sin(theta1 - theta2)))


# Calculate Q(V,theta)
def Q_calculado(v1, v2, theta1, theta2, g, b):
    return -v1*v1*b - (v1 * v2 * (g * math.sin(theta1 - theta2) - b * math.cos(theta1 - theta2)))


def read_dbar():
    #num;type;p;q;v;theta
    dbar_file = open('dbar1.txt')
    dbar_num = []
    dbar_type = []
    dbar_P = []
    dbar_Q = []
    dbar_V = []
    dbar_theta = []
    for line in dbar_file:
        splitLine = line.split(';')
        dbar_num.append(int(splitLine[0]))
        dbar_type.append(int(splitLine[1]))
        dbar_P.append(float(splitLine[2]))
        dbar_Q.append(float(splitLine[3]))
        dbar_V.append(float(splitLine[4]))
        dbar_theta.append(float(splitLine[5]))
    return dbar_num, dbar_type, dbar_P, dbar_Q, dbar_V, dbar_theta

def read_dlin():
    #de;para;r;x;bsh
    dlin_file = open('dlin1.txt')
    dlin_de = []
    dlin_para = []
    dlin_r = []
    dlin_x = []
    dlin_bsh = []
    for line in dlin_file:
        splitLine = line.split(';')
        dlin_de.append(int(splitLine[0]))
        dlin_para.append(int(splitLine[1]))
        dlin_r.append(float(splitLine[2]))
        dlin_x.append(float(splitLine[3]))
        dlin_bsh.append(float(splitLine[4]))
    return dlin_de, dlin_para, dlin_r, dlin_x, dlin_bsh


def get_rx_matrix(nbar, nlin):
    r_matrix = np.zeros((nbar, nbar))
    x_matrix = np.zeros((nbar, nbar))

    for i in range(nlin):
        de = dlin_de[i] - 1
        para = dlin_para[i] - 1
        r_matrix[de, para] = dlin_r[i]
        r_matrix[para, de] = r_matrix[de, para]

        x_matrix[de, para] = dlin_x[i]
        x_matrix[para, de] = x_matrix[de, para]

    return r_matrix, x_matrix

def get_gb_matrix(r_matrix, x_matrix):
    g_matrix = np.zeros((nbar, nbar))
    b_matrix = np.zeros((nbar, nbar))

    for i in range(nbar):
        for j in range(nbar):
            if i != j:
                g_matrix[i,j] = r_matrix[i, j]/(r_matrix[i, j] * r_matrix[i, j] + x_matrix[i, j]*x_matrix[i, j])
                b_matrix[i,j] = x_matrix[i, j]/(r_matrix[i, j]*r_matrix[i,j] + x_matrix[i, j]*x_matrix[i, j])

    return g_matrix, b_matrix

def get_bsh_matrix(nlin, dlin_de, dlin_para, dlin_bsh):
    bsh_matrix = np.zeros((nbar, nbar))
    for i in range(nlin):
        de = dlin_de[i] - 1
        para = dlin_para[i] - 1
        bsh_matrix[de,para] = dlin_bsh[i]/2
        bsh_matrix[para,de] = dlin_bsh[i]/2

    return bsh_matrix


def get_Y_matrix(g_matrix, b_matrix, bsh_matrix):
    G = np.zeros((nbar, nbar))
    B = np.zeros((nbar, nbar))
    for i in range(nbar):
        for j in range(nbar):
            if i != j:
                G[i,j] = -g_matrix[i,j]
                G[j,i] = G[i,j]
                B[i,j] = -b_matrix[i,j]
                B[j,i] = B[i,j]

    for i in range(nbar):
        somaG = 0
        somaB = 0
        for j in range(nbar):
            if (i != j):
                somaG += g_matrix[i, j]
                somaB += b_matrix[i, j] + bsh_matrix[i, j]
                    
        G[i, i] = somaG
        B[i, i] = somaB 

    return G, B
########################################################################
dbar_num, dbar_type, dbar_P, dbar_Q, dbar_V, dbar_theta = read_dbar()
dlin_de, dlin_para, dlin_r, dlin_x, dlin_bsh = read_dlin()

nbar = len(dbar_num)
nlin = len(dlin_de)

r_matrix, x_matrix = get_rx_matrix(nbar, nlin)

g_matrix, b_matrix = get_gb_matrix(r_matrix, x_matrix)

bsh_matrix = get_bsh_matrix(nlin, dlin_de, dlin_para, dlin_bsh)

G, B = get_Y_matrix(g_matrix, b_matrix, bsh_matrix)