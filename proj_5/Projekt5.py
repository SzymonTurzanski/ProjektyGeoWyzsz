from math import *

from numpy import *
from prettytable import PrettyTable

global kx0, ky0, kz0, kex, key, kez, kH

ag = 6378137
e2g = 0.00669437999013

ak = 6378245
e2k = 0.0066934215520398155
kx0 = -33.4297
ky0 = 146.5746
kz0 = 76.2865
kH = 0.8407728 * 0.000001
kex = (-0.35867 / 3600)
key = (-0.05283 / 3600)
kez = (0.84354 / 3600)

PhiA = radians(50.25)
LamA = radians(20.75)
PhiB = radians(50.00)
LamB = radians(20.75)
PhiC = radians(50.25)
LamC = radians(21.25)
PhiD = radians(50.00)
LamD = radians(21.25)
PhiS = radians(50 + 7/60 + 30.97362/3600)
LamS = radians(21 + 0/60 + 02.34392/3600)
PhiSS = radians(50 + 7/60 + 30/3600)
LamSS = radians(21 + 0/60 + 0/3600)



def toXYZ(fi, lambd, H, a, e2):
    N = a / sqrt(1 - e2 * (sin(fi) ** 2))
    x = (N + H) * cos(fi) * cos(lambd)
    y = (N + H) * cos(fi) * sin(lambd)
    z = (N * (1 - e2) + H) * sin(fi)
    return x, y, z


def conversion(deg_decimal):
    s = int(deg_decimal)
    min = int((deg_decimal - s) * 60)
    sec = (deg_decimal - s - min / 60) * 3600
    sec = "{:.5f}".format(sec)
    sec = str(sec)
    s, min = str(s), str(min)

    return s + '° ' + min + "' " + sec + "'' "


def hirvonen(x, y, z, a, e2):
    r = sqrt((x ** 2) + (y ** 2))
    fi1 = atan((z/r) * ((1-e2) ** -1))
    n = a / sqrt(1 - e2 * sin(fi1) ** 2)
    h = (r / cos(fi1)) - n
    fi2 = atan((z / r) * (1 - e2 * (n / (n + h))) ** -1)
    epsilon = radians(0.00005/3600)

    while abs(fi2-fi1) > epsilon:
        fi1 = fi2
        n = a / sqrt(1 - e2 * sin(fi1) ** 2)
        h = (r / cos(fi1)) - n
        fi2 = atan((z / r) * (1 - e2 * (n / (n + h))) ** -1)

    n = a / sqrt(1 - e2 * sin(fi2) ** 2)
    h = (r / cos(fi2)) - n
    lam = atan(y/x)

    fi_d = degrees(fi2)
    lam_d = degrees(lam)
    return conversion(fi_d), conversion(lam_d), h


def transform(xp, yp, zp):
    vector_m = array([[kx0], [ky0], [kz0]])
    coordinate_m = array([[xp], [yp], [zp]])
    rotation_m = array([[kH, radians(kez), radians(-key)], [radians(-kez), kH, radians(kex)], [radians(key), radians(-kex), kH]])
    m_0 = coordinate_m + rotation_m.dot(coordinate_m) + vector_m
    return transpose(m_0)

print("")
print("Współrzędne Punktu A w xyz w GRS80 ", (toXYZ(PhiA, LamA, 0, ag, e2g)))
print("")
print("Współrzędne Punktu B w xyz w GRS80", (toXYZ(PhiB, LamB, 0, ag, e2g)))
print("")
print("Współrzędne Punktu C w xyz w GRS80", (toXYZ(PhiC, LamC, 0, ag, e2g)))
print("")
print("Współrzędne Punktu D w xyz w GRS80", (toXYZ(PhiD, LamD, 0, ag, e2g)))
print("")
print("Współrzędne Punktu S w xyz w GRS80", (toXYZ(PhiS, LamS, 0, ag, e2g)))
print("")
print("Współrzędne Punktu SS w xyz w GRS80", (toXYZ(PhiSS, LamSS, 0, ag, e2g)))
print("")
print("Współrzędne Punktu A w xyz w e.Krasowskiego", (transform(3821451.635636481, 1447818.510760937, 4880617.05983242)))
print("")
print("Współrzędne Punktu B w xyz w e.Krasowskiego", (transform(3841408.3482922805, 1455379.43282719, 4862789.037706472)))
print("")
print("Współrzędne Punktu C w xyz w e.Krasowskiego", (transform(3808671.6868384206, 1481111.4156221119, 4880617.05983242)))
print("")
print("Współrzędne Punktu D w xyz w e.Krasowskiego", (transform(3828561.6589489407, 1488846.2027530423, 4862789.037706472)))
print("")
print("Współrzędne Punktu S w xyz w e.Krasowskiego", (transform(3825068.929880957, 1468306.3937127043, 4871714.592016324)))
print("")
print("Współrzędne Punktu SS w xyz w e.Krasowskiego", (transform(3825030.691234842, 1468341.5865585238, 4871733.878207174)))
print("")
print("Współrzędne Punktu A w e.Krasowskiego", (hirvonen(3821428.58996192, 1447942.18763557, 4880698.98862972, ak, e2k)))
print("")
print("Współrzędne Punktu B w e.Krasowskiego", (hirvonen(3841385.34575167, 1455503.06544474, 4862870.95955055, ak, e2k)))
print("")
print("Współrzędne Punktu C w e.Krasowskiego", (hirvonen(3808648.7665734, 1481235.17275336, 4880699.04979542, ak, e2k)))
print("")
print("Współrzędne Punktu D w e.Krasowskiego", (hirvonen(3828538.78247279, 1488969.91604633, 4862871.02103567, ak, e2k)))
print("")
print("Współrzędne Punktu S w e.Krasowskiego", (hirvonen(3825045.96875475, 1468430.08850005, 4871796.54802818, ak, e2k)))
print("")
print("Współrzędne Punktu SS w e.Krasowskiego", (hirvonen(3825007.73022534, 1468465.28149831, 4871815.83430624, ak, e2k)))


tabela = PrettyTable(['Punkty', 'Phi GRS80', 'Lambda GRS80', 'H GRS80', 'X GRS80', 'Y GRS80', 'Z GRS80', 'X Krasowski', 'Y Krasowski', 'Z Krasowski', 'Współrzędne Geodezyjne Krasowski Phi', 'Współrzędne Geodezyjne Krasowski Lambda', 'Współrzędne Geodezyjne Krasowski H'])
tabela.add_row(['A', "50°15’", "20°45'", '0', '3821451.636', '1447818.511', '4880617.060', '3821428.590', '1447942.188', '4880698.989', "50° 15' 1.05526''", "20° 45' 6.24968''", " -32.498"])
tabela.add_row(['B', "50°00'", "20°45'", '0', '3841408.348', '1455379.433', '4862789.038', '3841385.346', '1455503.065', '4862870.959', "50°00'01.06167''", "20°45'06.21437''", "-32.630"])
tabela.add_row(['C', "50°15'", "21°15'", '0', '3808671.687', '1481111.416', '4880617.060', '3808648.767', '1481235.173', '4880699.050', "50°15'01.02250''", "21°15'06.24111''", "-31.664"])
tabela.add_row(['D', '50°', "21°15'", '0', ' 3828561.659', ' 1488846.203', '4862789.038', '3828538.783', '1488969.916', '4862871.021 ', "50°00'01.03265''", "21°15'06.20584''", "-31.792"])
tabela.add_row(['S', "50°07'30.97362''", "21°00'02.34392''", '0', '33825030.691', '1468341.587', '4871733.878', '3825045.969', '1468430.089', '4871796.548', "50°07'32.01568''", "21°00'08.57170''", " -32.145"])
tabela.add_row(['SS', "50°07'30.0''", "21°00'00.0''", '0', '3825068.930', '1468306.394', '4871714.592', '3825007.730', '1468465.281', '4871815.834', " 50°07'31.04211", "21°00'06.22775''", " -32.146 "])

print(tabela)