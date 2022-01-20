from math import *

from numpy import *

# punkt sredniej szerokosci

phiSS = radians(50 + 7/60 + 30/3600)
lambdaSS = radians(21 + 0/60 + 0/3600)

# punkt srodkowy

phiS = radians(50 + 7/60 + 30.97362/3600)
lambdaS = radians(21 + 0/60 + 02.34392/3600)

# punkt A
phiA = radians(50.25)
lambdaA = radians(20.75)

# punkt B
phiB = radians(50)
lambdaB = radians(20.75)

# punkt C
phiC = radians(50.25)
lambdaC = radians(21.25)

# punkt D
phiD = radians(50)
lambdaD = radians(21.25)

a = 6378137
e2 = 0.00669437999013
b = a * sqrt(1 - e2)
eprim2 = (a **2 - b **2)/b **2


def liczN(fi):
    return a / sqrt(1 - e2 * (sin(fi)) ** 2)

def liczM(fi):
    return a * (1 - e2) / (sqrt(1 - e2 * (sin(fi)) ** 2)) ** 3

def GaussKruger(fi, la):
    L0 = 19 * (pi / 180)
    A0 = 1 - e2 / 4 - 3 * e2 ** 2 / 64 - 5 * e2 ** 3 / 256
    A2 = 3 * (e2 + e2 ** 2 / 4 + 15 * e2 ** 3 / 128) / 8
    A4 = 15 * (e2 ** 2 + 3 * e2 ** 3 / 4) / 256
    A6 = 35 * e2 ** 3 / 3072
    sigma = a * (A0 * fi - A2 * sin(2 * fi) + A4 * sin(4 * fi) - A6 * sin(6 * fi))
    t = tan(fi)
    n2 = eprim2 * cos(fi) ** 2
    l = la - L0
    Xgk = sigma + 0.5 * l ** 2 * liczN(fi) * sin(fi) * cos(fi) * (
                1 + l ** 2 / 12 * cos(fi) ** 2 * (5 - t ** 2 + 9 * n2 + 4 * n2 ** 2) + l ** 4 / 360 * cos(fi) ** 4 * (
                    61 - 58 * t ** 2 + t ** 4 + 270 * n2 - 330 * n2 * t ** 2))
    Ygk = l * liczN(fi) * cos(fi) * (1 + l ** 2 / 6 * cos(fi) ** 2 * (1 - t ** 2 + n2) + l ** 4 / 120 * cos(fi) ** 4 * (
                5 - 18 * t ** 2 + t ** 4 + 14 * n2 - 58 * n2 * t ** 2))
    Xgk = round(Xgk, 3)
    Ygk = round(Ygk, 3)

    return Xgk, Ygk

gkA = GaussKruger(phiA, lambdaA)
print("gkA" + str(gkA))

gkB = GaussKruger(phiB, lambdaB)
print("gkB" + str(gkB))

gkC = GaussKruger(phiC, lambdaC)
print("gkC" + str(gkC))

gkD = GaussKruger(phiD, lambdaD)
print("gkD" + str(gkD))

gkSS = GaussKruger(phiSS, lambdaSS)
print("gkSS" + str(gkSS))

gkS = GaussKruger(phiS, lambdaS)
print("gkS" + str(gkS))


def FLto2000(fi, lambd):
    if degrees(lambd) < 16.5:
        strefa = 15
    elif 16.5 <= degrees(lambd) < 19.5:
        strefa = 18
    elif 19.5 <= degrees(lambd) < 22.5:
        strefa = 21
    elif degrees(lambd) >= 22.5:
        strefa = 24

    L0 = strefa * (pi/180)
    N = a / (sqrt(1 - e2 * (sin(fi)) ** 2))
    A0 = 1 - (e2 / 4) - ((3 * (e2) ** 2) / 64) - ((5 * (e2) ** 3) / 256)
    A2 = (3 / 8) * (e2 + ((e2) ** 2) / 4 + ((15 * (e2) ** 3) / 128))
    A4 = (15 / 256) * ((e2) ** 2 + (3 * (e2) ** 3) / 4)
    A6 = ((35 * (e2) ** 3) / 3072)
    t = tan(fi)
    L = lambd - L0
    n2 = eprim2 * (cos(fi) ** 2)
    sigma = a * (A0 * fi - A2 * sin(2 * fi) + A4 * sin(4 * fi) - A6 * sin(6 * fi))

    x_GK = sigma + (L ** 2 / 2) * N * sin(fi) * cos(fi) * (1 + (L ** 2 / 12) * (cos(fi) ** 2) * (5 - (t ** 2) + 9 * n2 + 4 * (n2 ** 2)) + ((L ** 4) / 360) * (cos(fi) ** 4) * (61 - 58 * (t ** 2) + (t ** 4) + 270 * n2 - 330 * n2 * (t ** 2)))
    y_GK = L * N * cos(fi) * (1 + ((L ** 2) / 6) * (cos(fi) ** 2) * (1 - (t ** 2) + n2) + ((L ** 4) / 120) * (cos(fi) ** 4) * (5 - 18 * (t ** 2) + (t ** 4) + 14 * n2 - 58 * n2 * (t ** 2)))

    m0 = 0.999923
    nr_strefy = strefa/3
    x_2000 = m0 * x_GK
    y_2000 = m0 * y_GK + 1000000 * nr_strefy + 500000
    x_2000 = round(x_2000, 3)
    y_2000 = round(y_2000, 3)
    return x_2000, y_2000,

FLto2000A = FLto2000 (phiA, lambdaA)
print("2000A" + str(FLto2000A))

FLto2000B = FLto2000 (phiB, lambdaB)
print("2000B" + str(FLto2000B))

FLto2000C = FLto2000 (phiC, lambdaC)
print("2000C" + str(FLto2000C))

FLto2000D = FLto2000 (phiD, lambdaD)
print("2000D" + str(FLto2000D))

FLto2000SS = FLto2000 (phiSS, lambdaSS)
print("2000SS" + str(FLto2000SS))

FLto2000S = FLto2000 (phiS, lambdaS)
print("2000S" + str(FLto2000S))




def FLto92(fi, la):
    Xgk, Ygk = GaussKruger(fi, la)
    x92 = 0.9993 * Xgk - 5300000
    y92 = 0.9993 * Ygk + 500000
    x92 = round(x92, 3)
    y92 = round(y92, 3)

    return x92, y92



FLto1992A = FLto92 (phiA, lambdaA)
print("1992A" + str(FLto1992A))

FLto1992B = FLto92 (phiB, lambdaB)
print("1992B" + str(FLto1992B))

FLto1992C = FLto92 (phiC, lambdaC)
print("1992C" + str(FLto1992C))

FLto1992D = FLto92 (phiD, lambdaD)
print("1992D" + str(FLto1992D))

FLto1992SS = FLto92 (phiSS, lambdaSS)
print("1992SS" + str(FLto1992SS))

FLto1992S = FLto92 (phiS, lambdaS)
print("1992S" + str(FLto1992S))



def angles_formatting(degrees):
    deg = int(degrees)
    minutes = int((degrees - deg)*60)
    seconds = (degrees - deg - minutes/60) * 3600
    seconds = round(seconds, 5)
    deg = f"{deg:02d}"
    minutes = f"{minutes:02d}"
    int_sec = f"{int(seconds):02d}"
    float_seconds = str(round(seconds-int(seconds), 5))[1:]

    return f'{deg}°{minutes}\'{int_sec}{float_seconds}"'



def GKtoFL(x, y):
    A0 = 1 - e2/4 - 3*e2**2/64 - 5*e2**3/256
    A2 = 3*(e2 + e2**2/4 + 15*e2**3/128)/8
    A4 = 15*(e2**2 + 3*e2**3/4)/256
    A6 = 35*e2**3/3072

    mianow = a*A0
    fi0 = x/mianow

    while True:
        sigma = a * (A0 * fi0 - A2 * sin(2 * fi0) + A4 * sin(4 * fi0) - A6 * sin(6 * fi0))
        fi1 = fi0 + (x - sigma) / a * A0
        if abs(fi1 - fi0) < radians(0.000001 / 3600):
            break
        else:
            fi0 = fi1

    N = liczN(fi1)
    M = liczM(fi1)
    t = tan(fi1)
    n2 = eprim2 * cos(fi1) ** 2
    L0 = 21

    fi = fi1 - y**2*t/(2*M*N)*(1 - y**2/(12*N**2)*(5+3*t**2+n2-9*n2*t**2-4*n2**2) + y**4/(360*N**4)*(61+90*t**2+45*t**4))
    la = L0*pi/180 + y/(N*cos(fi1))*(1 - y**2/(6*N**2)*(1+2*t**2+n2) + y**4/(120*N**4)*(5+28*t**2+24*t**4+6*n2+8*n2*t**2))
    # fi = angles_formatting(fi)
    # la = angles_formatting(la)

    return fi, la






def u92toFL(x, y):
    Xgk = (x+5300000)/0.9993
    Ygk = (y-500000)/0.9993
    x, y = GKtoFL(Xgk, Ygk)
    return x, y


def u2ktoFL(x, y):
    Xgk = x / 0.999923
    L0 = int((y-500000)/1000000)*3
    Ygk = (y - 500000) % 1000000 / 0.999923
    x, y = GKtoFL(Xgk, Ygk)
    return x, y

def skala_1992(x, y):
    a = 6378137
    e2 = 0.0066943800290
    fi, lam = u92toFL(x, y)
    M = a * (1 - e2)/(1 - e2 * sin(fi) ** 2) ** (3/2)
    N = a / (1 - e2 * sin(fi) ** 2) ** 0.5

    Q = sqrt(M * N)

    mgk = 1 + y ** 2/(2 * Q ** 2) + y ** 2/(24 * Q ** 4)
    m92 = 0.9993 * mgk
    kappa = (1 - m92)*1000

    return m92, kappa


def skala_2000(x, y):
    m0 = 0.999923
    fi, lam = u2ktoFL(x, y)
    M = a * (1 - e2) / (1 - e2 * sin(fi) ** 2) ** (3 / 2)
    N = a / (1 - e2 * sin(fi) ** 2) ** 0.5
    Q = sqrt(M * N)

    mgk = 1 + y ** 2 / (2 * Q ** 2) + y ** 2 / (24 * Q ** 4)
    m2000 = m0 * mgk
    kappa = (1 - m2000)*1000

    return m2000, kappa


def GK_skala(ygk, phi):
    # phi, lamb = GKtoFL(xgk, ygk)
    a = 6378137
    e2 = 0.00669437999013
    M = (a * (1 - e2)) / (sqrt((1 - e2 * (sin(phi) ** 2)) ** 3))
    N = a / (sqrt(1 - e2 * (sin(phi) ** 2)))
    R = sqrt(M * N)
    m = 1 + (ygk ** 2) / (2 * R ** 2) + (ygk ** 4) / (24 * R ** 4)
    kappa = (m - 1) * 1000
    return m, kappa

odwGKA = GKtoFL(gkA[0], gkA[1])
odwGKB = GKtoFL(gkB[0], gkB[1])
odwGKC = GKtoFL(gkC[0], gkC[1])
odwGKD = GKtoFL(gkD[0], gkD[1])
odwGKSS = GKtoFL(gkSS[0], gkSS[1])
odwGKS = GKtoFL(gkS[0], gkS[1])

odw92A = u92toFL(FLto1992A[0], FLto1992A[1])
odw92B = u92toFL(FLto1992B[0], FLto1992B[1])
odw92C = u92toFL(FLto1992C[0], FLto1992C[1])
odw92D = u92toFL(FLto1992D[0], FLto1992D[1])
odw92SS = u92toFL(FLto1992SS[0], FLto1992SS[1])
odw92S = u92toFL(FLto1992S[0], FLto1992S[1])

odw2000A = u2ktoFL(FLto2000A[0], FLto2000A[1])
odw2000B = u2ktoFL(FLto2000B[0], FLto2000B[1])
odw2000C = u2ktoFL(FLto2000C[0], FLto2000C[1])
odw2000D = u2ktoFL(FLto2000D[0], FLto2000D[1])
odw2000SS = u2ktoFL(FLto2000SS[0], FLto2000SS[1])
odw2000S = u2ktoFL(FLto2000S[0], FLto2000S[1])




GKA_skala = GK_skala(gkA[1], odwGKA[0])
print("GKA_skala", GKA_skala)

GKB_skala = GK_skala(gkB[1], odwGKB[0])
print("GKB_skala", GKB_skala)

GKC_skala = GK_skala(gkC[1], odwGKC[0])
print("GKC_skala", GKC_skala)

GKD_skala = GK_skala(gkD[1], odwGKD[0])
print("GKD_skala", GKD_skala)

GKSS_skala = GK_skala(gkSS[1], odwGKSS[0])
print("GKSS_skala", GKSS_skala)

GKS_skala = GK_skala(gkS[1], odwGKS[0])
print("GKS_skala", GKS_skala)

skala_1992A = GKA_skala[0] * 0.9993
skala_1992A = [skala_1992A, (skala_1992A - 1) * 1000]
print("skala_1992A", skala_1992A)

skala_1992B = GKB_skala[0] * 0.9993
skala_1992B = [skala_1992B, (skala_1992B - 1) * 1000]
print("skala_1992B", skala_1992B)

skala_1992C = GKC_skala[0] * 0.9993
skala_1992C = [skala_1992C, (skala_1992C - 1) * 1000]
print("skala_1992C", skala_1992C)

skala_1992D = GKD_skala[0] * 0.9993
skala_1992D = [skala_1992D, (skala_1992D - 1) * 1000]
print("skala_1992D", skala_1992D)

skala_1992SS = GKSS_skala[0] * 0.9993
skala_1992SS = [skala_1992SS, (skala_1992SS - 1) * 1000]
print("skala_1992SS", skala_1992SS)

skala_1992S = GKS_skala[0] * 0.9993
skala_1992S = [skala_1992S, (skala_1992S - 1) * 1000]
print("skala_1992S", skala_1992S)


skala_2000A = GKA_skala[0] * 0.999923
skala_2000A = [skala_2000A, (skala_2000A - 1) * 1000]
print("skala_2000A", skala_2000A)

skala_2000B = GKB_skala[0] * 0.999923
skala_2000B = [skala_2000B, (skala_2000B - 1) * 1000]
print("skala_2000B", skala_2000B)

skala_2000C = GKC_skala[0] * 0.999923
skala_2000C = [skala_2000C, (skala_2000C - 1) * 1000]
print("skala_2000C", skala_2000C)

skala_2000D = GKD_skala[0] * 0.999923
skala_2000D = [skala_2000D, (skala_2000D - 1) * 1000]
print("skala_2000D", skala_2000D)

skala_2000SS = GKSS_skala[0] * 0.999923
skala_2000SS = [skala_2000SS, (skala_2000SS - 1) * 1000]
print("skala_2000SS", skala_2000SS)

skala_2000S = GKS_skala[0] * 0.999923
skala_2000S = [skala_2000S, (skala_2000S - 1) * 1000]
print("skala_2000S", skala_2000S)


m2_GKA = [GKA_skala[0] ** 2, (GKA_skala[0] ** 2 - 1) * 10000]
print ("m2_GKA", m2_GKA)
m2_GKB = [GKB_skala[0] ** 2, (GKB_skala[0] ** 2 - 1) * 10000]
print ("m2_GKB", m2_GKB)
m2_GKC = [GKC_skala[0] ** 2, (GKC_skala[0] ** 2 - 1) * 10000]
print ("m2_GKC", m2_GKC)
m2_GKD = [GKD_skala[0] ** 2, (GKD_skala[0] ** 2 - 1) * 10000]
print ("m2_GKD", m2_GKD)
m2_GKSS = [GKSS_skala[0] ** 2, (GKSS_skala[0] ** 2 - 1) * 10000]
print ("m2_GKSS", m2_GKSS)
m2_GKS = [GKS_skala[0] ** 2, (GKS_skala[0] ** 2 - 1) * 10000]
print ("m2_GKS", m2_GKS)

m2_skala_2000A = [skala_2000A[0] ** 2, (skala_2000A[0] ** 2 - 1) * 10000]
print ("m2_skala_2000A", m2_skala_2000A)
m2_skala_2000B = [skala_2000B[0] ** 2, (skala_2000B[0] ** 2 - 1) * 10000]
print ("m2_skala_2000B", m2_skala_2000B)
m2_skala_2000C = [skala_2000C[0] ** 2, (skala_2000C[0] ** 2 - 1) * 10000]
print ("m2_skala_2000C", m2_skala_2000C)
m2_skala_2000D = [skala_2000D[0] ** 2, (skala_2000D[0] ** 2 - 1) * 10000]
print ("m2_skala_2000D", m2_skala_2000D)
m2_skala_2000SS = [skala_2000SS[0] ** 2, (skala_2000SS[0] ** 2 - 1) * 10000]
print ("m2_skala_2000SS", m2_skala_2000SS)
m2_skala_2000S = [skala_2000S[0] ** 2, (skala_2000S[0] ** 2 - 1) * 10000]
print ("m2_skala_2000S", m2_skala_2000S)

m2_skala_1992A = [skala_1992A[0] ** 2, (skala_1992A[0] ** 2 -1) * 10000]
print ("m2_skala_1992A", m2_skala_1992A)
m2_skala_1992B = [skala_1992B[0] ** 2, (skala_1992B[0] ** 2 -1) * 10000]
print ("m2_skala_1992B", m2_skala_1992B)
m2_skala_1992C = [skala_1992C[0] ** 2, (skala_1992C[0] ** 2 -1) * 10000]
print ("m2_skala_1992C", m2_skala_1992C)
m2_skala_1992D = [skala_1992D[0] ** 2, (skala_1992D[0] ** 2 -1) * 10000]
print ("m2_skala_1992D", m2_skala_1992D)
m2_skala_1992SS = [skala_1992SS[0] ** 2, (skala_1992SS[0] ** 2 -1) * 10000]
print ("m2_skala_1992SS", m2_skala_1992SS)
m2_skala_1992S = [skala_1992S[0] ** 2, (skala_1992SS[0] ** 2 -1) * 10000]
print ("m2_skala_1992S", m2_skala_1992S)



print(f" Pole elipsoidalne z poprzedniego zadania 994265196.080189")

print(f"pole Gaussa-Krugera ", abs((gkA[0] - gkD[0]) * (gkA[1] - gkD[1])))

print(f"pole układu 2000 ", abs((FLto2000A[0] - FLto2000D[0]) * (FLto2000A[1] - FLto2000D[1])))

print(f"pole układu 1992 ", abs((FLto1992A[0] - FLto1992D[0]) * (FLto1992A[1] - FLto1992D[1])))













