import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
import time
import datetime
import os
from shutil import copy
import main_gui
from tkinter import *
from ttest import *
from PyQt5.QtCore import QThread, pyqtSignal


class Reservoir():
    def __init__(self, Ps, rw, Rk, Nr, Nz):
        self.P_0 = 101325
        self.T_0 = 273.15
           #transform in (Pa) unit
        self.Ps = Ps*101325
        self.rw = rw
        self.Rk = Rk
        self.Nr = Nr
        self.Nz = Nz
        #self.m = m
        #self.k = k * 1E-15          #transform in (m2) unit
        #self.rho = rho
        #self.Cp = Cp
        #self.lyam = lyam
        #self.s1_0 = s1_0
        #[[0] * m for i in range(n)]
        self.P = [[0.0] * self.Nr for i in range(self.Nz)]
        self.P_pr = [[0.0] * self.Nr for i in range(self.Nz)]            # np.zeros(Nr)
        self.P_it = [[0.0] * self.Nr for i in range(self.Nz)]
        self.T = [[0.0] * self.Nr for i in range(self.Nz)]
        self.T_pr = [[0.0] * self.Nr for i in range(self.Nz)]
        self.T_it = [[0.0] * self.Nr for i in range(self.Nz)]
        self.r = [0.0] * self.Nr
        self.s1 = [[0.0] * self.Nr for i in range(self.Nz)]
        self.s1_pr = [[0.0] * self.Nr for i in range(self.Nz)]
        self.s1_it = [[0.0] * self.Nr for i in range(self.Nz)]
        self.k1 = [[0.0] * self.Nr for i in range(self.Nz)]
        self.k2 = [[0.0] * self.Nr for i in range(self.Nz)]
        self.Pres = [0.0] * self.Nz


class Wellbore():
    def __init__(self, Ph0, Ph1, dt, dz, Nz, Tb, Gz, zt, zb, s1_0, v1b, v2b, Tin, sigma):
        self.Ph0 = Ph0 * 101325
        self.Ph1 = Ph1 * 101325
        self.dt = dt
        self.dz = dz
        self.Nz = Nz
        self.Tb = Tb
        self.Gz = Gz
        self.zt = zt
        self.zb = zb
        self.s1_0 = s1_0
        self.v1b = (-1)*v1b
        self.v2b = (-1)*v2b
        self.Tin = Tin
        self.sigma = sigma
        self.z = [0.0] * self.Nz
        self.v1 = [0.0] * self.Nz
        self.v1_pr = [0.0] * self.Nz
        self.v1_it = [0.0] * self.Nz
        self.v2 = [0.0] * self.Nz
        self.v2_pr = [0.0] * self.Nz
        self.v2_it = [0.0] * self.Nz
        self.vm = [0.0] * self.Nz
        self.S1 = [0.0] * self.Nz
        self.S1_pr = [0.0] * self.Nz
        self.S1_it = [0.0] * self.Nz
        self.S2 = [0.0] * self.Nz
        self.S2_pr = [0.0] * self.Nz
        self.P = [0.0] * self.Nz
        self.P_pr = [0.0] * self.Nz
        self.P_it = [0.0] * self.Nz
        self.T = [0.0] * self.Nz
        self.T_pr = [0.0] * self.Nz
        self.T_it = [0.0] * self.Nz
        self.Tgeo = [0.0] * self.Nz
        self.rho1 = [0.0] * self.Nz
        self.rho1_pr = [0.0] * self.Nz
        self.rho2 = [0.0] * self.Nz
        self.rho2_pr = [0.0] * self.Nz
        self.C0 = 1.2
        self.Fw = 0.0
        self.S = 0.0
        self.it = 0
        self.dP = 0.0
        self.dS = 0.0
        self.dv1 = 0.0
        self.dv2 = 0.0
        self.dT = 0.0


class Fluid():
    def __init__(self, rho_0, Cp, lyam, mu, eps, nu, R, L12, beta):
        self.rho_0 = rho_0
        #self.rho = np.zeros(Reservoir.Nr)
        self.Cp = Cp
        self.lyam = lyam
        self.eps = eps/101325   #transform in (K/Pa) unit
        self.nu = nu/101325
        self.mu = mu * 1E-3     #transform in (Pa*s) unit
        self.R = R/101325      #объемный коэффициент растворимости
        self.L12 = L12 * 1E3    #transform in (J/kg) unit
        self.beta = beta/101325



class solver():
    def __init__(self, dt, time, dt_exp):
        self.dt = dt
        self.time = time
        self.dt_exp = dt_exp
        self.Nt = round(self.time // self.dt)


def load(file_name):

    f = open(file_name, 'r')
    f.readline()
    data = f.readline().split()       #Time properties
    global solv
    solv = solver(float(data[0]), int(data[1]), int(data[2]))

    f.readline()
    data1 = f.readline().split()        # Well radius and Reservoir
    f.readline()
    data2 = f.readline().split()        # Head pressure,atm Bubble point pressure,atm	  sigma, N/m
    f.readline()
    data3 = f.readline().split()        # Well top
    Nz = round((float(data3[1])-float(data3[0])) // float(data3[2])) + 1
    f.readline()
    data4 = f.readline().split()        # bottomhole temp
    f.readline()
    data5 = f.readline().split()        # wellbore init. oil holdup

    f.readline()
    data = f.readline().split()  # Oil properties
    global oil
    oil = Fluid(float(data[0]),float(data[1]),float(data[2]),float(data[3]),float(data[4]),float(data[5]),float(data[6]),0, float(data[7]))


    f.readline()

    data = f.readline().split()  # Gas properties
    global gas
    gas = Fluid(float(data[0]),float(data[1]),float(data[2]),float(data[3]),float(data[4]),float(data[5]),0,float(data[6]),0)

    global water
    water = Fluid(1100, 4170, 0.6, 1, 0.02, 0.0015, 0, 0, 0)

    f.readline()        # Well loggers depth, m
    global loggs
    loggs = f.readline().split()

    f.readline()        # Inflow intervals
    f.readline()        # Reservoir properties

    data = [line for line in f.readlines()]

    global plast
    plast = Reservoir(float(data2[3]), float(data1[0]), float(data1[1]), int(data1[2]), int(Nz))

    plast.ht = [0.0] * len(data)
    plast.hb = [0.0] * len(data)
    plast.k = [0.0] * len(data)
    plast.m = [0.0] * len(data)
    plast.rho = [0.0] * len(data)
    plast.Cp = [0.0] * len(data)
    plast.lyam = [0.0] * len(data)
    plast.S1_0 = [0.0] * len(data)

    global count_plast
    count_plast = len(data)

    oil.rho = [[0.0] * plast.Nr for i in range(Nz)]
    oil.rho_pr = [[0.0] * plast.Nr for i in range(Nz)]
    gas.rho = [[0.0] * plast.Nr for i in range(Nz)]

    global rho1, rho1_pr, rho2, rho2_pr, cgw, cgw_pr

    rho1 = [0.0] * Nz
    rho1_pr = [0.0] * Nz
    rho2 = [0.0] * Nz
    rho2_pr = [0.0] * Nz
    cgw = [0.0] * Nz
    cgw_pr = [0.0] * Nz

    global well
    well = Wellbore(float(data2[0]),float(data2[1]),float(data2[2]), float(data3[2]), Nz, float(data4[0]), float(data4[1]),
                    float(data3[0]), float(data3[1]), float(data5[0]), float(data5[1]), float(data5[2]), float(data5[3]), float(data2[4]))

    well.Fw = 2*math.pi*plast.rw
    well.S = math.pi * (plast.rw)**2

    if (len(data)==0):
        return 0

    k = data[0].split()

    dat = [[0.0] * len(k)] * len(data)

    for i in range(len(data)):
        dat[i] = data[i].split()

    plast.ht = get_column(0, dat, len(data))
    plast.hb = get_column(1, dat, len(data))
    plast.k = get_column(2, dat, len(data))

    for i in range(len(plast.k)):
        plast.k[i] = float(plast.k[i]) * 1E-15

    plast.m = get_column(3, dat, len(data))
    plast.rho = get_column(4, dat, len(data))
    plast.Cp = get_column(5, dat, len(data))
    plast.lyam = get_column(6, dat, len(data))
    plast.S1_0 = get_column(7, dat, len(data))

    f.close()


def get_column(col_cout, dat, fl_c):
    p = []
    for i in range(fl_c):
        p += [dat[i][col_cout]]
    return p


def sr(a,b):
    if (a!=b):
        return (a - b) / math.log(a/b)


def C_g(p):

#  defined by Henry's law

    if (p<plast.Ps):
        return 1/(1 + oil.rho_0/(gas.rho_0*oil.R*p))

    if (p>=plast.Ps):
        return 1/(1 + oil.rho_0/(gas.rho_0*oil.R*plast.Ps))


def OFP():
    for i in range(well.Nz):
        if (z_ind[i] == 1):
            for j in range(plast.Nr):
                if (plast.s1_pr[i][j]>0.15):
                    plast.k1[i][j] = ((plast.s1_pr[i][j]-0.15)/0.85 ) ** 2.8 * (3.4 - 2.4 * plast.s1_pr[i][j])
                else: plast.k1[i][j] = 0
                plast.k2[i][j] = ((1-plast.s1_pr[i][j]) ** 3.5) * (1 + 3 * plast.s1_pr[i][j])


def plast_def(z):
    i = 0
    res = 0.0
    for i in range(len(plast.ht)):
        if (z >= float(plast.ht[i]) ) and (z <= float(plast.hb[i])):
            res = 1.0
            return res
        else:
            res = 0.0
    return res

def zumpf(z):
    res = 0.0
    if (z > float(plast.hb[-1])):
        res = 1.0
        return res
    else:
        res = 0.0
    return res

def initial():
    well.z[0] = well.zt
    well.z[well.Nz-1] = well.zb
    well.T[well.Nz-1] = well.Tb
    well.Tgeo[well.Nz-1] = well.Tb
    well.T_pr[well.Nz-1] = well.Tb
    well.P[0] = well.Ph0
    plast.Pres[0] = well.Ph0

    global grav
    grav = 9.81

    global z_ind
    z_ind = [0.0] * plast.Nz
    global sump
    sump = [0.0] * well.Nz

    for i in range(well.Nz-1):
        well.z[i+1] = well.z[i] + well.dz
        z_ind[i+1] = plast_def(well.z[i+1])
        sump[i + 1] = zumpf(well.z[i + 1])
        if (sump[i + 1] != 1):
            well.P[i + 1] = well.P[0] + oil.rho_0 * grav * (well.z[i + 1] - well.z[0])
        else:
            well.P[i + 1] = well.P[i] + water.rho_0 * grav * well.dz

        plast.Pres[i+1] = well.P[i+1]

    for i in range(well.Nz-2, -1, -1 ):
        well.T[i] = well.T[i+1] - well.Gz * well.dz
        well.T_pr[i] = well.T[i]
        well.Tgeo[i] = well.T[i]

    for i in range(well.Nz):
        well.P_pr[i] = well.P[i]
        well.P_it[i] = well.P[i]
        well.T_it[i] = well.T[i]
        well.S1[i] = well.s1_0
        well.S1_pr[i] = well.s1_0
        well.S2[i] = 1 - well.s1_0
        well.S2_pr[i] = 1 - well.s1_0


    plast.r[0] = plast.rw
    plast.r12 = [0.0] * plast.Nr
    plast.dr2 = [0.0] * plast.Nr

    oil.c_g = [[0.0] * plast.Nr for i in range(well.Nz)]
    oil.c_g_pr = [[0.0] * plast.Nr for i in range(well.Nz)]


    for i in range(plast.Nr-1):
        plast.r[i+1] = plast.r[i] * (plast.Rk/plast.rw) ** (1/(plast.Nr - 1))
        plast.r12[i] = sr(plast.r[i + 1], plast.r[i])
        if i > 0: plast.dr2[i] = (plast.r12[i] ** 2 - plast.r12[i - 1] ** 2) / 2

    for i in range(well.Nz):
        if (z_ind[i]==1):
            for j in range(plast.Nr):
                plast.T_pr[i][j] = well.T[i]
                plast.T_it[i][j] = well.T[i]
                plast.T[i][j] = well.T[i]
                plast.P[i][j] = well.P[i]
                plast.P_pr[i][j] = well.P[i]
                plast.P_it[i][j] = well.P[i]
                plast.s1[i][j] = well.S1[i]
                plast.s1_it[i][j] = well.S1[i]
                plast.s1_pr[i][j] = well.S1[i]
                oil.c_g[i][j] = C_g(plast.P[i][j])
                oil.c_g_pr[i][j] = oil.c_g[i][j]

    well.v1 = [0.0] * plast.Nz
    well.v1_pr = [0.0] * plast.Nz
    well.v2 = [0.0] * plast.Nz
    well.v2_pr = [0.0] * plast.Nz

    global J1, J2

    J1 = [0.0] * plast.Nz
    J2 = [0.0] * plast.Nz

    oil.u = [[0.0] * plast.Nr for i in range(well.Nz)]
    gas.u = [[0.0] * plast.Nr for i in range(well.Nz)]

    global K_r, m_p, rho_r, lyam_r, Cp_r, S1_in_r
    K_r = [0.0] * well.Nz
    m_p = [0.0] * well.Nz
    rho_r = [0.0] * well.Nz
    lyam_r = [0.0] * well.Nz
    Cp_r = [0.0] * well.Nz
    S1_in_r = [0.0] * well.Nz

    i = 0

    for m in range(len(plast.ht)):
        for i in range(well.Nz):
            if (well.z[i] >= float(plast.ht[m])) and (well.z[i] <= float(plast.hb[m])):
                K_r[i] = plast.k[m]
                m_p[i] = float(plast.m[m])
                rho_r[i] = float(plast.rho[m])
                lyam_r[i] = float(plast.lyam[m])
                Cp_r[i] = float(plast.Cp[m])
                S1_in_r[i] = float(plast.S1_0[m])

    global log_j
    log_j = [0] * (len(loggs)+1)
    i=1
    for j in range(well.Nz):
        if (i>(len(loggs))):
            break
        if (float(loggs[i-1])==well.z[j]):
            log_j[i] = j
            i += 1


    OFP()


def next_iter():
    for i in range(well.Nz):
        if (z_ind[i] == 1):
            for j in range(plast.Nr):
                plast.P_it[i][j] = plast.P[i][j]
                plast.T_it[i][j] = plast.T[i][j]
                plast.s1_it[i][j] = plast.s1[i][j]
                oil.c_g[i][j] = C_g(plast.P[i][j])


def next_step():
    for i in range(well.Nz):
        well.T_pr[i] = well.T[i]
        well.T_it[i] = well.T[i]
        well.P_pr[i] = well.P[i]
        well.P_it[i] = well.P[i]
        well.S1_pr[i] = well.S1[i]
        well.S1_it[i] = well.S1[i]
        well.S2_pr[i] = well.S2[i]
        well.v1_pr[i] = well.v1[i]
        well.v1_it[i] = well.v1[i]
        well.v2_pr[i] = well.v2[i]
        well.v2_it[i] = well.v2[i]

        if (z_ind[i] == 1):
            for j in range(plast.Nr):
                plast.P_pr[i][j] = plast.P[i][j]
                plast.P_it[i][j] = plast.P[i][j]
                plast.T_pr[i][j] = plast.T[i][j]
                plast.T_it[i][j] = plast.T[i][j]
                plast.s1_pr[i][j] = plast.s1[i][j]
                plast.s1_it[i][j] = plast.s1[i][j]
                oil.c_g_pr[i][j] = oil.c_g[i][j]
    OFP()


def res_dens_fl():
    for i in range(well.Nz):
        if (z_ind[i] == 1):
            for j in range(plast.Nr):
                oil.rho[i][j] = oil.rho_0 * (1 + oil.beta*(plast.P[i][j] - plast.P_0))
                oil.rho_pr[i][j] = oil.rho_0 * (1 + oil.beta*(plast.P_pr[i][j] - plast.P_0))
                gas.rho[i][j] = gas.rho_0 * plast.P[i][j] / plast.P_0 #* 300 / plast.T_pr[i][j]


def well_dens():
    for j in range(well.Nz):
        rho1[j] = oil.rho_0 * (1 + oil.beta * (well.P[j] - plast.P_0))
        rho1_pr[j] = oil.rho_0 * (1 + oil.beta * (well.P_pr[j] - plast.P_0))
        rho2[j] = gas.rho_0 * well.P[j] / plast.P_0 #* plast.T_0 / well.T[j]
        rho2_pr[j] = gas.rho_0 * well.P_pr[j] / plast.P_0 #* plast.T_0 / well.T_pr[j]


def progonka(a,b,c,d):
    n = len(a)
    alpha = [0.0]*n
    beta = [0.0]*n
    x = [0.0]*n
    i=0
    for i in range(n-1):
        if (b[i]!=0):
            alpha[i+1] = - c[i] / (b[i] + a[i] * alpha[i])
            beta[i+1] = (-a[i]*beta[i] + d[i]) / (b[i] + a[i] * alpha[i])

    x[n-1] = d[n-1]

    for i in reversed(range(n-1)):
        x[i] = alpha[i+1] * x[i+1] + beta[i+1]

    return x

def res_pressure_solve(j):
    a = [0.0]*plast.Nr  #np.zeros(plast.Nr)
    b = [0.0]*plast.Nr
    c = [0.0]*plast.Nr
    d = [0.0]*plast.Nr

    for i in range(1, plast.Nr - 1):
        a[i] = (1 - oil.c_g[j][i])/gas.rho[j][i] * plast.r12[i-1]*K_r[j]*plast.k1[j][i]/oil.mu * 1/(plast.r[i] - plast.r[i-1])
        a[i] = a[i] + (1-oil.c_g[j][i])/oil.rho[j][i]*plast.r12[i-1]*K_r[j]*plast.k2[j][i]/gas.mu * 1/(plast.r[i] - plast.r[i-1])
        a[i] = a[i] - (1/gas.rho[j][i] - 1/oil.rho[j][i])*plast.r12[i-1]*(1-oil.c_g[j][i])*K_r[j]*plast.k1[j][i]/oil.mu * 1/(plast.r[i] - plast.r[i-1])

        b[i] = -(1-oil.c_g[j][i])/gas.rho[j][i]*plast.r12[i]*K_r[j]*plast.k1[j][i+1]/oil.mu * 1/(plast.r[i+1] - plast.r[i])
        b[i] = b[i] - (1-oil.c_g[j][i])/gas.rho[j][i]*plast.r12[i-1]*K_r[j]*plast.k1[j][i]/oil.mu * 1/(plast.r[i] - plast.r[i-1])
        b[i] = b[i] - (1-oil.c_g[j][i])/oil.rho[j][i]*plast.r12[i]*K_r[j]*plast.k2[j][i+1]/gas.mu * 1/(plast.r[i+1] - plast.r[i])
        b[i] = b[i] - (1-oil.c_g[j][i])/oil.rho[j][i]*plast.r12[i-1]*K_r[j]*plast.k2[j][i]/gas.mu * 1/(plast.r[i] - plast.r[i-1])
        b[i] = b[i] + (1/gas.rho[j][i] - 1/oil.rho[j][i])*plast.r12[i] * (1-oil.c_g[j][i+1])*K_r[j]*plast.k1[j][i+1]/oil.mu * 1/(plast.r[i+1] - plast.r[i])
        b[i] = b[i] + (1/gas.rho[j][i] - 1/oil.rho[j][i])*plast.r12[i-1] * (1-oil.c_g[j][i])*K_r[j]*plast.k1[j][i]/oil.mu * 1/(plast.r[i] - plast.r[i-1])
        b[i] = b[i] - plast.dr2[i]/solv.dt*m_p[j]*(1-oil.c_g[j][i])/gas.rho[j][i]*(plast.s1[j][i]*oil.beta + (1-plast.s1[j][i])/oil.rho[j][i]*gas.rho_0/plast.P_0)
        b[i] = b[i] + plast.dr2[i]/solv.dt*m_p[j]*plast.s1[j][i]*oil.beta*(1-oil.c_g[j][i])*(1/gas.rho[j][i] - 1/oil.rho[j][i])

        c[i] = (1-oil.c_g[j][i])/gas.rho[j][i]*plast.r12[i]*K_r[j]*plast.k1[j][i+1]/oil.mu * 1/(plast.r[i+1] - plast.r[i])
        c[i] = c[i] + (1-oil.c_g[j][i])/oil.rho[j][i]*plast.r12[i]*K_r[j]*plast.k2[j][i+1]/gas.mu * 1/(plast.r[i+1] - plast.r[i])
        c[i] = c[i] - (1/gas.rho[j][i] - 1/oil.rho[j][i])*plast.r12[i] * (1-oil.c_g[j][i+1])*K_r[j]*plast.k1[j][i+1]/oil.mu * 1/(plast.r[i+1] - plast.r[i])

        d[i] = - plast.dr2[i]/solv.dt*m_p[j]*(1-oil.c_g[j][i])/gas.rho[j][i]*(plast.s1[j][i]*oil.beta + (1-plast.s1[j][i])/oil.rho[j][i]*gas.rho_0/plast.P_0)*plast.P_pr[j][i]
        d[i] = d[i] + plast.dr2[i]/solv.dt*m_p[j]*plast.s1[j][i]*oil.beta*(1-oil.c_g[j][i])*(1/gas.rho[j][i] - 1/oil.rho[j][i]) * plast.P_pr[j][i]

    a[0] = 0
    b[0] = 1
    c[0] = 0
    d[0] = well.P[j]

    a[plast.Nr - 1] = 0
    b[plast.Nr - 1] = 1
    c[plast.Nr - 1] = 0
    d[plast.Nr - 1] = plast.Pres[j]

    plast.P[j] = progonka(a,b,c,d)


def gas_concetr():
    for j in range(well.Nz):
        if (z_ind[j] == 1):
            for i in range(plast.Nr):
                oil.c_g[j][i] = C_g(plast.P_it[j][i])
                oil.c_g_pr[j][i] = C_g(plast.P_pr[j][i])


def saturation_solve(j):
    plast.s1[j][plast.Nr - 1] = S1_in_r[j]

    for i in reversed(range(1, plast.Nr - 1)):
        plast.s1[j][i] = (1 / ((1 - oil.c_g[j][i])*oil.rho[j][i])) * ((1 - oil.c_g_pr[j][i]) * oil.rho_pr[j][i] * plast.s1_pr[j][i]  + solv.dt / (m_p[j] * plast.dr2[i]) * oil.rho[j][i] *
                                                (plast.r12[i] * (1 - oil.c_g[j][i + 1]) * K_r[j] * plast.k1[j][i + 1] / oil.mu * (plast.P[j][i + 1] - plast.P[j][i]) / (plast.r[i + 1] - plast.r[i]) -
                                                plast.r12[i - 1] * (1 - oil.c_g[j][i]) * K_r[j] * plast.k1[j][i] / oil.mu * (plast.P[j][i] - plast.P[j][i - 1]) / (plast.r[i] - plast.r[i - 1])))

        if (plast.s1[j][i]>1):
            plast.s1[j][i] = 1
        elif (plast.s1[j][i]<0):
            plast.s1[j][i] = 0

    plast.s1[j][0] = plast.s1[j][1]


def reservoir_phase_vel(j):
    for i in range(plast.Nr - 1):
        oil.u[j][i] = - K_r[j]*plast.k1[j][i+1]/oil.mu * (plast.P[j][i+1] - plast.P[j][i])/(plast.r[i+1] - plast.r[i])
        gas.u[j][i] = - K_r[j]*plast.k2[j][i+1]/gas.mu * (plast.P[j][i+1] - plast.P[j][i])/(plast.r[i+1] - plast.r[i])


def psr(a,b):
    if (a==0) and (b==0):
        return 0
    else:
        return  2 * a * b / (a + b)


def lyam_eff(j,i):
    lyam_fl = (oil.lyam ** (m_p[j]*plast.s1[j][i])) * (gas.lyam ** (m_p[j]*(1 - plast.s1[j][i])))
    return (lyam_fl**m_p[j]) * (lyam_r[j]**(1-m_p[j]))


def C_res(j,i):
    return m_p[j] * (oil.rho[j][i]*plast.s1[j][i]*oil.Cp + gas.rho[j][i]*(1-plast.s1[j][i])*gas.Cp) + (1 - m_p[j])*rho_r[j]*Cp_r[j]


def J_12r(j,i):
    if (plast.P[j][i]>=plast.Ps):
        return 0
    else:
        return -m_p[j] * (oil.rho[j][i]*plast.s1[j][i] - oil.rho_pr[j][i]*plast.s1_pr[j][i]) / solv.dt - \
               oil.rho[j][i]/(plast.r[i]*(plast.r12[i]-plast.r12[i-1])) * (plast.r12[i]*oil.u[j][i] - plast.r12[i-1]*oil.u[j][i-1])


def J_12w(j):
    if (sump[j]==1):
        return 0
    elif (well.P[j]<plast.Ps):
        return J1[j] - (well.S1[j]*rho1[j] - well.S1_pr[j]*rho1_pr[j])/solv.dt - ( - well.S1[j]*rho1[j]*well.v1[j] + well.S1[j+1]*rho1[j+1]*well.v1[j+1])/well.dz
    elif (well.P[j]>=plast.Ps):
        return 0


def q_DT(j,i):
    #return - (oil.eps*oil.Cp*oil.rho[j][i])*(plast.r12[i]*oil.u[j][i]*(plast.P[j][i+1]-plast.P[j][i])/(plast.r[i+1]-plast.r[i]) - plast.r12[i-1]*oil.u[j][i-1]*(plast.P[j][i]-plast.P[j][i-1])/(plast.r[i]-plast.r[i-1]))
    return -1/2*(oil.eps*oil.rho[j][i]*oil.Cp*K_r[j]*plast.k1[j][i]/oil.mu + gas.eps*gas.Cp*gas.rho[j][i]*K_r[j]*plast.k2[j][i]/gas.mu) * \
          ( ((plast.P[j][i]-plast.P[j][i-1])/(plast.r[i]-plast.r[i-1]))**2*(plast.r[i]**2 - plast.r12[i-1]**2) +
            ((plast.P[j][i+1]-plast.P[j][i])/(plast.r[i+1]-plast.r[i]))**2*(plast.r12[i]**2 - plast.r[i] ** 2) )


def q_AD(j,i):
    return -m_p[j]*(plast.r12[i]**2-plast.r12[i-1]**2)/2*(oil.rho[j][i]*oil.Cp*plast.s1[j][i]*oil.nu + gas.rho[j][i]*gas.Cp*(1-plast.s1[j][i])*gas.nu)*(plast.P[j][i]-plast.P_pr[j][i])


def res_calc_temp(j):
    a = [0.0]*plast.Nr
    b = [0.0]*plast.Nr
    c = [0.0]*plast.Nr
    d = [0.0]*plast.Nr

    for i in range(1, plast.Nr - 1):
        l1 = lyam_eff(j,i-1)
        l2 = lyam_eff(j,i)
        l3 = lyam_eff(j,i+1)
        a[i] = solv.dt *psr(l2,l1)*2*plast.r12[i-1]/(plast.r[i]-plast.r[i-1])
        b[i] = -(solv.dt*psr(l2,l1) * 2*plast.r12[i-1]/(plast.r[i]-plast.r[i-1]) +
                 solv.dt*psr(l3,l2)*2*plast.r12[i]/(plast.r[i+1]-plast.r[i]) + plast.dr2[i]*C_res(j,i) -
                 solv.dt*plast.r[i]*(oil.rho[j][i+1]*oil.Cp*oil.u[j][i] + gas.rho[j][i+1]*gas.Cp*gas.u[j][i]))
        c[i] = solv.dt*psr(lyam_eff(j,i+1),lyam_eff(j,i))*2*plast.r12[i]/(plast.r[i+1]-plast.r[i]) - solv.dt*plast.r[i]*(oil.rho[j][i+1]*oil.Cp*oil.u[j][i] +
                                                                                                                     gas.rho[j][i+1]*gas.Cp*gas.u[j][i])
        d[i] = -plast.dr2[i]*C_res(j,i)*plast.T_pr[j][i] - plast.dr2[i]*solv.dt*J_12r(j,i)*gas.L12 + q_DT(j,i)*solv.dt + q_AD(j,i)

    a[0] = 0
    b[0] = 1
    c[0] = -1
    d[0] = 0

    a[plast.Nr - 1] = 0
    b[plast.Nr - 1] = 1
    c[plast.Nr - 1] = 0
    d[plast.Nr - 1] = well.Tgeo[j]

    plast.T[j] = progonka(a, b, c, d)


def Re(j):
    return 2*plast.rw * well.vm[j] * (well.S1[j]*rho1[j] + (1-well.S1[j])*rho2[j])/(well.S1[j]*oil.mu + (1-well.S1[j])*gas.mu)


def tau(j):
    f = 0.0
    if (0<Re(j)<2300):
        f = 16/Re(j)
    elif (2300<Re(j)<1E5):
        f = 0.079 * (Re(j))**(-0.25)
    elif (Re(j)==0.0):
        return 0.0
    return 0.5 * f * (well.S1[j]*rho1[j] + (1-well.S1[j])*rho2[j]) * well.vm[j] * math.fabs(well.vm[j])


def well_press_solve():
    #    well.P[j + 1] = well.P[j] + (well.S1[j]*rho1[j] + (1-well.S1[j])*rho2[j]) * grav * well.dz #+ 2*tau(j)/plast.rw * well.dz
    for j in range(well.Nz - 1):
        if (sump[j] != 1):
            well.P[j + 1] = well.P[j] + (well.S1[j] * rho1[j] + (1 - well.S1[j]) * rho2[j]) * grav * well.dz + tau(j) * well.Fw / well.S * well.dz
        else:
            well.P[j + 1] = well.P[j] + water.rho_0 * grav * well.dz


def well_holdup():
    well.S1[well.Nz-1] = well.s1_0
    well.S2[well.Nz-1] = 1 - well.s1_0

    for j in reversed(range(1, well.Nz - 1)):
        if (sump[j] != 1):
            well.S1[j] = 1 / (well.dz/solv.dt*rho1[j]*(1-C_g(well.P[j])) - well.v1[j]*rho1[j]*(1-C_g(well.P[j]))) * \
                         (well.dz*J1[j]*(1-C_g(well.P[j])) + well.dz/solv.dt*rho1_pr[j]*(1-C_g(well.P_pr[j]))*well.S1_pr[j] -
                          well.v1[j+1]*rho1[j+1]*(1-C_g(well.P[j+1]))*well.S1[j+1])
        if (well.S1[j]>1.0): well.S1[j]=1.0
        if (well.S1[j] < 0.0): well.S1[j] = 0.0
        well.S2[j] = 1 - well.S1[j]

    well.S1[0]=well.S1[1]
    well.S2[0]=well.S2[1]


def v_slip(j):
    u_c = (well.sigma * grav*(rho1[j]-rho2[j])/(rho1[j]**2))**0.25
    return -1.53*u_c


def well_velocities():
    well.v1[well.Nz-1] = well.v1b
    well.v2[well.Nz-1] = well.v2b

    for j in reversed(range(1, well.Nz - 1)):
        if (sump[j] != 1):
            well.v1[j] = 1/((-1)*well.S1[j]*rho1[j])* (-well.S1[j+1]*rho1[j+1]*well.v1[j+1] + (J1[j]+J2[j])*well.dz -
                                                   (well.S2[j+1]*rho2[j+1]*well.v2[j+1]-well.S2[j]*rho2[j]*well.v2[j]) -
                                                   well.dz/solv.dt*( (well.S1[j]*rho1[j]+well.S2[j]*rho2[j]) - (well.S1_pr[j]*rho1_pr[j] + well.S2_pr[j]*rho2_pr[j]) ))

            well.vm[j] = well.S1[j]*well.v1[j] + well.S2[j]*well.v2[j]
            well.v2[j] = well.C0*well.vm[j] + v_slip(j)

    well.v1[0] = well.v1[1]
    well.vm[0] = well.vm[1]
    well.v2[0] = well.v2[1]


def zn_v1(j):
    if (well.v1[j]<=0):
        return 1
    elif (well.v1[j]>0):
        return 0


def lyam_m(j):
    return well.S1[j]*oil.lyam + well.S2[j]*gas.lyam


def Nusselt(j):
    prdl = 0.0
    Re1 = 2100
    Re2 = 4000

    if (oil.lyam!=0) and (gas.lyam!=0):
        prdl = (well.S1[j]*oil.Cp + well.S2[j]*gas.Cp)*(well.S1[j]*oil.mu + well.S2[j]*gas.mu) / lyam_m(j)
    else:   prdl = 0.0

    Rnld = 2 * plast.rw * math.fabs(well.vm[j]) * (well.S1[j]*rho1[j]+well.S2[j]*rho2[j]) / (well.S1[j]*oil.mu + well.S2[j]*gas.mu)

    if (sump[j] == 1):
        prdl = water.Cp*water.mu/water.lyam
        Rnld = 0.0

    if (Rnld<=Re1):
        return 4.36
    elif (Re1<Rnld<Re2):
        return 4.36+(0.021*Re2**0.8*prdl**0.43 - 4.36)*(Rnld-Re1)/(Re2-Re1)
    elif (Rnld>=4000):
        return 0.021*Rnld**0.8*prdl**0.43


def alpha(j):
    if (sump[j] != 1):
        return lyam_m(j)*Nusselt(j) / (2*plast.rw)
    else:
        return water.lyam * Nusselt(j) / (2*plast.rw)

def H_obm(j,ttime):
    rcas_i = 0.1
    rcas_o = 0.105
    rcem_i = 0.105
    rcem_o = 0.125
    lyam_cas = 40
    lyam_cem = 1.3
    lyam_rr = 1.5
    ro_r = 2500
    c_r = 500

    t_d = lyam_rr * ttime / (ro_r * c_r * rcem_o * rcem_o)
    if (t_d <= 1.5):
        f_t = 1.1281 * math.sqrt(t_d) * (1 - 0.3 * math.sqrt(t_d))
    else:   f_t = (0.4063 + 0.5 * math.log(t_d)) * (1 + 0.6 / t_d)

    return 2 / (plast.rw ) * ((rcas_i * (1 / (rcas_i * alpha(j)) + math.log(rcas_o / rcas_i) / lyam_cas) + math.log(rcem_o / rcem_i) / lyam_cem + f_t / lyam_rr) ** (-1))


def well_temp_calc(tt):
    a = [0.0] * well.Nz
    b = [0.0] * well.Nz
    c = [0.0] * well.Nz
    d = [0.0] * well.Nz

    for j in range(1, well.Nz-1):
        a[j] = -lyam_m(j)/(well.dz**2)
        #a[j] += -(1-zn_v1(j+1))* (well.S1[j]*oil.Cp*rho1[j]*well.v1[j+1] + well.S2[j]*gas.Cp*rho2[j]*well.v2[j+1])

        b[j] = well.dz/solv.dt*(well.S1[j]*oil.Cp*rho1[j] + well.S2[j]*gas.Cp*rho2[j])
        b[j] += -zn_v1(j+1)*(well.S1[j+1]*oil.Cp*rho1[j+1]*well.v1[j+1] + well.S2[j+1]*gas.Cp*rho2[j+1]*well.v2[j+1])
        b[j] += (1 - zn_v1(j+1)) * (well.S1[j] * oil.Cp * rho1[j] * well.v1[j+1] + well.S2[j]*gas.Cp*rho2[j]*well.v2[j+1] )
        b[j] += H_obm(j,tt)*well.dz*(1-z_ind[j]) + well.dz*( oil.Cp*J1[j] + gas.Cp*J2[j] ) + 2 * lyam_m(j) / (well.dz ** 2)

        c[j] = zn_v1(j+1)*(well.S1[j+1]*oil.Cp*rho1[j+1]*well.v1[j+1] + well.S2[j+1]*gas.Cp*rho2[j+1]*well.v2[j+1]) - lyam_m(j)/(well.dz**2)
        c[j] += -(1 - zn_v1(j+1)) * (well.S1[j] * oil.Cp * rho1[j] * well.v1[j+1] + well.S2[j]*gas.Cp*rho2[j]*well.v2[j+1] )

        d[j] = H_obm(j,tt)*well.Tgeo[j]*well.dz*(1-z_ind[j]) + well.dz/solv.dt*(well.S1[j]*oil.Cp*rho1[j] + well.S2[j]*gas.Cp*rho2[j])*well.T_pr[j] + \
               (oil.Cp*J1[j] + gas.Cp*J2[j])*well.dz*plast.T[j][0] + gas.L12*math.fabs(J_12w(j))*well.dz*(1-z_ind[j])

        if (sump[j] != 1):
            d[j] += (oil.Cp*rho1[j]*well.S1[j]*oil.nu+gas.Cp*rho2[j]*well.S2[j]*gas.nu)*(well.P[j]-well.P_pr[j])*well.dz/solv.dt
        else:
            d[j] += water.rho_0*water.Cp*water.nu*(well.P[j]-well.P_pr[j])*well.dz/solv.dt

    a[0] = 0
    b[0] = 1
    c[0] = - 1
    d[0] = 0

    a[well.Nz-1] = 0
    b[well.Nz-1] = 1
    c[well.Nz-1] = 0
    d[well.Nz-1] = well.Tin

    well.T = progonka(a,b,c,d)


def jacobian(f, x):
    n = 5
    Jac = zeros([n, n])
    f0 = f(x)
    for i in arange(0, n, 1):
        tt = x[i]
        x[i] = tt + h
        f1 = f(x)
        x[i] = tt
        Jac[:, i] = (f1 - f0) / h
    return Jac, f0


def f(x):
    f = zeros([n])
    for i in arange(0, n - 1, 1):
        f[i] = (5 + 2 * math.sin(x[i])) * x[i] - x[i - 1] - 2 * x[i + 1] - 2
    f[0] = (3 + 2 * x[0]) * x[0] - 2 * x[1] - 3
    f[n - 1] = (3 + 2 * x[n - 1]) * x[n - 1] - x[n - 2] - 4
    return f


#def newton_method():


def get_contin():

    for j in range(well.Nz):
        if (z_ind[j] == 1):
            plast.dP[j] = math.fabs(plast.P[j][0]-plast.P_it[j][0])
            plast.dS[j] = math.fabs(plast.s1[j][0]-plast.s1_it[j][0])
            plast.dT[j] = math.fabs(plast.T[j][0]-plast.T_it[j][0])

            for i in range(plast.Nr):
                if (plast.dP[j]<math.fabs(plast.P[j][i]-plast.P_it[j][i])):
                    plast.dP[j] = math.fabs(plast.P[j][i]-plast.P_it[j][i])
                if (plast.dS[j]<math.fabs(plast.s1[j][i]-plast.s1_it[j][i])):
                    plast.dS[j] = math.fabs(plast.s1[j][i]-plast.s1_it[j][i])
                if (plast.dT[j]<math.fabs(plast.T[j][i]-plast.T_it[j][i])):
                    plast.dT[j] = math.fabs(plast.T[j][i]-plast.T_it[j][i])


def mass_sources():
    for j in range(well.Nz):
        if (z_ind[j] == 1):
            J1[j] = 2 / plast.rw * oil.rho[j][0]*math.fabs(oil.u[j][0])
            J2[j] = 2 / plast.rw * gas.rho[j][0] * math.fabs(gas.u[j][0])


def set_well_press(t):
    k1 = 0.005
    tyst = 120

    well.P[0] = well.Ph1 + (well.Ph0 - well.Ph1) * math.exp(-k1*t)

    #if (t <= tyst):
    #    well.P[0] = well.Ph0 - (well.Ph0 - well.Ph1) * t / tyst
    #else:
    #    well.P[0] = well.Ph1


def reservoir_solve():
    eps_P = 1E-5
    eps_S = 1E-5
    eps_T = 1E-8

    plast.dP = [0.0] * well.Nz
    plast.dS = [0.0] * well.Nz
    plast.dT = [0.0] * well.Nz

    plast.dP[0] = well.P[0]
    plast.dS[0] = well.S1[0]
    plast.dT[0] = well.T[0]
    plast.it = 0
    max_iter = 200

    res_dens_fl()

    while (max(plast.dP)>eps_P) or (max(plast.dS)>eps_S) or (max(plast.dT)>eps_T):
        plast.dP[0] = 0.0
        plast.dS[0] = 0.0
        plast.dT[0] = 0.0

        for j in range(well.Nz):
            if (z_ind[j] == 1):
                res_pressure_solve(j)
                saturation_solve(j)
                reservoir_phase_vel(j)
                res_calc_temp(j)


        gas_concetr()
        #graph_r(100)
        res_dens_fl()
        get_contin()
        next_iter()
        plast.it += 1
        if plast.it>max_iter:
            print('iteration count above max number!')
            exit()
        #print(max(plast.dP))
    mass_sources()
    print ('Reservoir task was solved. Iters number = ' + str(plast.it))
    inf =  'Reservoir task was solved. Iters number = ' + str(plast.it)
    return inf


def get_cont_w():
    for j in range(well.Nz):
        well.dP = math.fabs(well.P[j] - well.P_it[j])
        well.dS = math.fabs(well.S1[j] - well.S1_it[j])
        well.dv1 = math.fabs(well.v1[j] - well.v1_it[j])
        well.dv2 = math.fabs(well.v2[j] - well.v2_it[j])
        well.dT = math.fabs(well.T[j] - well.T_it[j])


def next_iter_w():
    for j in range(well.Nz):
        well.P_it[j] = well.P[j]
        well.S1_it[j] = well.S1[j]
        well.T_it[j] = well.T[j]
        well.v1_it[j] = well.v1[j]
        well.v2_it[j] = well.v2[j]


def wellbore_solve(t):
    eps_P = 1E-6
    eps_S = 1E-5
    eps_v1 = 1E-5
    eps_v2 = 1E-5
    eps_T = 1E-8

    max_iter = 500
    well_dens()

    well.dP = well.P[0]
    well.dS = well.S1[0]
    well.dT = well.T[0]
    well.dv1 = 1
    well.dv2 = 1

    well.it=0

    while (well.dP > eps_P) or (well.dS > eps_S) or (well.dT > eps_T) or (well.dv1 > eps_v1) or (well.dv2 > eps_v2):
        well_press_solve()
        well_holdup()
        well_velocities()
        well_temp_calc(t)
        well.it += 1
        get_cont_w()
        next_iter_w()
        well_dens()
        #print('dP='+str(well.dP)+' dS='+str(well.dS)+' dv1='+str(well.dv1)+' dv2='+str(well.dv2))
        if well.it>max_iter:
            print('iteration count above max number!')
            exit()
    print('Wellbore task was solved. Iters number = ' + str(well.it))
    inf = 'Wellbore task was solved. Iters number = ' + str(well.it)
    return inf


def graph(t, j):
    fig = plt.figure()
    p_ax = fig.add_subplot(3,1,1)
    p_ax.set_title('Pressure')
    s_ax = fig.add_subplot(3,1,2)
    s_ax.set_title('Saturation')
    t_ax = fig.add_subplot(3,1,3)
    t_ax.set_title('Temperature')
    p_ax.plot(plast.r,plast.P[j])
    s_ax.plot(plast.r,plast.s1[j])
    t_ax.plot(plast.r,plast.T[j])
    plt.savefig(str(t)+'.png')
    #plt.show()


def graph_w_v1(t,catal):
    plt.gca().invert_yaxis()
    plt.title('Oil velocity     time = ' + str(t)+'s')
    plt.xlabel('velocity, m/s')
    plt.ylabel('depth, m')
    plt.grid(color='gray', linestyle='dashed')
    plast_cout = len(plast.ht)
    for j in range(plast_cout):
        plt.axhspan(float(plast.ht[j]), float(plast.hb[j]), facecolor='0.5', alpha=0.5)
    plt.plot(well.v1, well.z, label='oil')
    plt.savefig(catal+'/'+'Oil-Vel-well_' + str(t) + '.png')
    plt.draw()
    plt.close()
    #plt.show()


def graph_w_v2(t,catal):
    plt.gca().invert_yaxis()
    plt.title('Gas velocity     time = ' + str(t)+'s')
    plt.xlabel('velocity, m/s')
    plt.ylabel('depth, m')
    plt.grid(color='gray', linestyle='dashed')
    plast_cout = len(plast.ht)
    for j in range(plast_cout):
        plt.axhspan(float(plast.ht[j]), float(plast.hb[j]), facecolor='0.5', alpha=0.5)
    plt.plot(well.v2, well.z, label='gas')
    plt.savefig(catal+'/'+'Gas-Vel-well_' + str(t) + '.png')
    plt.draw()
    plt.close()


def graph_w_T(t,catal):
    plt.gca().invert_yaxis()
    plt.title('Well temperature,K   time = ' + str(t)+'s')
    plt.xlabel('T, K')
    plt.ylabel('depth, m')
    plt.grid(color='gray', linestyle='dashed')
    plast_cout = len(plast.ht)
    for j in range(plast_cout):
        plt.axhspan(float(plast.ht[j]), float(plast.hb[j]), facecolor='0.5', alpha=0.5)
    plt.plot(well.Tgeo, well.z, label='geotherm')
    plt.plot(well.T, well.z, label='well temper')
    plt.legend()
    plt.savefig(catal+'/'+'Temp-well_' + str(t) + '.png')
    plt.draw()
    plt.close()
    #plt.show()


def graph_w_S2(t,catal):
    plt.gca().invert_yaxis()
    plt.title('Well gas volume fraction,   time = ' + str(t)+'s')
    plt.xlabel('S2, ')
    plt.ylabel('depth, m')
    plt.grid(color='gray', linestyle='dashed')
    plast_cout = len(plast.ht)
    for j in range(plast_cout):
        plt.axhspan(float(plast.ht[j]), float(plast.hb[j]), facecolor='0.5', alpha=0.5)
    plt.plot(well.S2, well.z)
    plt.savefig(catal+'/'+'S2-well_' + str(t) + '.png')
    plt.draw()
    plt.close()
    #plt.show()


def graph_r(j,t,catal):
    plt.plot(plast.r, plast.T[j])
    plt.title('Reservoir temperature, z='+ str(well.z[j]) +  'm, time = ' + str(t)+'s')
    plt.xlabel('Radius, m')
    plt.ylabel('Temperature, K')
    plt.grid(color='gray', linestyle='dashed')
    plt.savefig(catal+'/'+'Temp-res_' + str(well.z[j]) +'m_'+ str(t) + 's.png')
    plt.draw()
    plt.close()
    #plt.show()


def export_r(t,catal):
    jj = 0
    if len(plast.ht)!= 0:
        for j in range(well.Nz-1):
            if (math.fabs(well.z[j] - float(plast.ht[0])) < 1E-5):
                jj = j
    else:
        return 0
    file = catal +'/'+ 'res_'+str(well.z[jj])+'m_'+str(t)+'s.txt'
    outfile = open(file, 'w')
    outfile.write('# r,m\tP,atm\ts1\tu1\tu2\tT,grad.C\n')
    x = plast.r
    y = plast.P[jj]
    a = plast.s1[jj]
    b = oil.u[jj]
    c = gas.u[jj]
    d = plast.T[jj]
    for xi, yi, ai, bi, ci, di in zip(x, y, a, b, c, d):
        outfile.write('%10.5f\t%10.5f\t%10.5f\t%10.7f\t%10.5f\t%10.7f\n' % (float(xi),float(yi)/101325,float(ai),float(bi),float(ci),float(di)))
    outfile.close()


def export_w(t,catal):
    file = catal +'/'+ 'well_'+str(t)+'s.txt'
    outfile = open(file, 'w')
    outfile.write('# z,m\tP,atm\ts1\tv1\tv2\tT,grad.C\n')
    x = well.z
    y = well.P
    a = well.S1
    b = well.v1
    c = well.v2
    d = well.T
    for xi, yi, ai, bi, ci, di in zip(x, y, a, b, c, d):
        outfile.write('%10.5f\t%10.5f\t%10.12f\t%10.10f\t%10.10f\t%10.5f\n' % (float(xi),float(yi),float(ai),float(bi),float(ci),float(di)))
    outfile.close()


def logger_well_point(time,catal,j):
    if time == solv.dt:
        outfile = open(catal + '/'+'point-well-'+str(well.z[j])+'m.txt', 'w')
        outfile.write('# time,hr\tP,Pa\tT,K\ts2\tv1,m/s\tv2,m/s\tQ1,m3/d\tQ2,m3/d\tres_iters\twell_iters\n')
        outfile.close()

    Q1 = math.pi * plast.rw **2 * math.fabs(well.v1[j]) * well.S1[j] * 86400
    Q2 = math.pi * plast.rw **2 * math.fabs(well.v2[j]) * well.S2[j] * 86400
    P = well.P[j]
    T = well.T[j]
    S2 = well.S2[j]
    u1 = well.v1[j]
    u2 = well.v2[j]
    it_r = plast.it
    it_w = well.it
    outfile = open(catal + '/' + 'point-well-' + str(well.z[j]) + 'm.txt', 'a')
    outfile.write('%10.5f\t%10.6f\t%10.6f\t%10.12f\t%10.10f\t%10.10f\t%10.5f\t%10.5f' % ( (time/3600), P, T, S2, u1, u2, Q1, Q2))
    outfile.write('\t' + str(it_w) + '\n')  # + str(it_r) + '\t'
    outfile.close()


def graph_surf():

    z = np.linspace(well.zt, well.zb, well.Nz)#np.meshgrid(well.z)
    r = plast.r
    #np.linspace(plast.rw, plast.Rk, plast.Nr) #np.meshgrid(plast.r)

    dat = np.zeros((well.Nz,plast.Nr))

    min1 = 0.0
    min2 = 1E6
    max1 = 0.0
    max2 = 0.0

    for i in range(well.Nz):
        if (z_ind[i] == 1):
            for j in range(plast.Nr):
                dat[i,j] = plast.P[i][j]
            min1 = min(plast.P[i])
            if (min2>min1): min2 = min1
            max1 = max(plast.P[i])
            if (max2<max1): max2 = max1

    v = np.linspace(min2, max2, 20, endpoint=True)

    plt.gca().invert_yaxis()
    #plt.contour(r,z,dat,colors='k')
    plt.contourf(r, z, dat, v, cmap=plt.cm.jet)
    plt.ylim(float(plast.hb[0]), float(plast.ht[0]))
    plt.colorbar(ticks=v)
    #plt.show()


def loggers(tt, cat):
    for i in range(len(log_j)):
        logger_well_point(tt,cat,log_j[i])

def stop_calc():
    main_gui.stop = True

def calculate(win):
    start_time = time.time()
    global k
    k = 1
    param_file = 'param.txt'
    load(param_file)
    initial()
    ttime = 0.0
    #dt_gr = 60

    win.stop = False
    now = datetime.datetime.now()
    dir_name = now.strftime("%m-%d-%Y %H-%M-%S")  # %H:%M
    path0=os.getcwd()
    path1=path0.replace('\\', '/')

    win.countbar = External()
    #win.countbar.countChanged.connect(win.onCountChanged)
    win.countbar.start()

    newpath =path1+'/'+dir_name
    log_files_path = newpath + '/' + 'loggers files'
    graph_path = newpath + '/' + 'plots'
    if not os.path.exists(newpath):
        os.mkdir(newpath)
    if not os.path.exists(log_files_path):
        os.mkdir(log_files_path)
    if not os.path.exists(graph_path):
        os.mkdir(graph_path)

    copy(path1+'/'+param_file, newpath+'/'+param_file)

    while (ttime < solv.time):
        ttime += solv.dt
        dt_inf = int(win.edittext1.text())
        set_well_press(ttime)
        info=''
        info = reservoir_solve()
        info += '\n'
        info += ' ' + wellbore_solve(ttime)
        info += '\n'
        loggers(ttime, newpath)
        next_step()
        if (round(ttime) % solv.dt_exp) < 1E-5 :
            export_w(ttime, log_files_path)
            export_r(ttime, log_files_path)

        print('Time = ' + str(ttime) + ' s')
        info += ' ' + 'Time = ' + str(ttime) + ' s' +'\n'

        win.text.appendPlainText(info)
        val = round(ttime * 100 / solv.time)
        #win.countbar.countChanged.connect()
        #win.onCountChanged(val)

        if (round(ttime % dt_inf)) < 1E-5 :
            win.figure.clear()
            win.plot(well.P, well.z)
            win.plot1(well.S1, well.z)
            win.plot2(well.T, well.z)
            win.plot_vel(well.v1, well.v2, well.z)
            win.plot_well(win.ht, win.hb)


        if win.stop == True:
            break


    info = ''
    print('Calculation complete.')
    print("-- %s seconds --" % (time.time() - start_time))
    info += 'Calculation complete.\n'
    info += "-- %s seconds --" % (time.time() - start_time)
    win.text.appendPlainText(info)



if __name__ == '__main__':
    b = 0





