from math import *

from numpy import *

s = 51123
nr = 10
e2 = 0.00669438002290
a = 6378137
b = a * sqrt(1 - e2)

# wspolrzedne punktow
phiA = (50.25)
lambdaA = (20.75)


phiD = (50)
lambdaD = (21.25)

# dane wierzchołków
A = [phiA, lambdaA]
D = [phiD, lambdaD]

phiSrednie = (phiA + phiD)/2
lamdaSrednie = (lambdaA + lambdaD)/2
Mean = [phiSrednie, lamdaSrednie]


def vincent(A, B):
    A = [radians(i) for i in A]
    B = [radians(i) for i in B]

    f = 1 - (b/a)
    delta_lambda = B[1] - A[1]
    L = delta_lambda
    Ua = arctan((1 - f)*tan(A[0]))
    Ub = arctan((1 - f)*tan(B[0]))

    while True:
        sin_sigma = sqrt((cos(Ub) * sin(L)) ** 2 + (cos(Ua) * sin(Ub) - sin(Ua) * cos(Ub) * cos(L)) ** 2)
        cos_sigma = sin(Ua) * sin(Ub) + cos(Ua) * cos(Ub) * cos(L)
        sigma = arctan(sin_sigma / cos_sigma)

        sin_a = (cos(Ua)*cos(Ub)*sin(L))/sin_sigma
        cos2_a = 1 - (sin_a**2)
        cos2_sigma_m = cos_sigma - (2*sin(Ua)*sin(Ub))/cos2_a
        C = (f / 16)*cos2_a*(4 + f*(4-3*cos2_a))
        L_i = delta_lambda + (1 - C)*f*sin_a*(sigma + C*sin_sigma*(cos2_sigma_m+C*cos_sigma*(1-2*(cos2_sigma_m**2))))

        if fabs(radians(L_i - L)) < (0.000001 / 3600):
            break
        else:
            L = L_i

    u2 = (a**2 - b**2)*cos2_a/(b**2)
    A = 1 + (u2/16384) * (4096 + u2*(-768+u2*(320 - 175*u2)))
    B = (u2/1024) * (256 + u2*(-128 + u2*(74 - 47*u2)))

    delta_sigma = B*sin_sigma*(cos2_sigma_m + 0.25*B*(cos_sigma*(-1+2*(cos2_sigma_m**2)) - 1/6*B*cos2_sigma_m*(-3 + 4*(sin_sigma**2))*(-3+4*(cos2_sigma_m**2))))

    s_AB = b*A*(sigma - delta_sigma)

    # obliczanie azymutu
    X_A_ab = cos(Ub)*sin(L)
    Y_A_ab = cos(Ua)*sin(Ub) - sin(Ua)*cos(Ub)*cos(L)
    A_ab = arctan(X_A_ab/Y_A_ab)
    if X_A_ab > 0 and Y_A_ab > 0:
        pass
    elif X_A_ab > 0 and Y_A_ab < 0:
        A_ab = A_ab + pi
    elif X_A_ab < 0 and Y_A_ab < 0:
        A_ab = A_ab + pi
    elif X_A_ab < 0 and Y_A_ab > 0:
        A_ab = A_ab + 2*pi

    X_A_ba = cos(Ua)*sin(L)
    Y_A_ba = -sin(Ua)*cos(Ub) + cos(Ua)*sin(Ub)*cos(L)
    A_ba = arctan(X_A_ba/Y_A_ba)
    if X_A_ba > 0 and Y_A_ba > 0:
        pass
    elif X_A_ba > 0 and Y_A_ba < 0:
        A_ba = A_ba + 2*pi
    elif X_A_ba < 0 and Y_A_ba < 0:
        A_ba = A_ba + 2*pi
    elif X_A_ba < 0 and Y_A_ba > 0:
        A_ba = A_ba + 3*pi
    return s_AB, degrees(A_ab), degrees(A_ba)


def kivioji():
    s_AB, A_ab = vincent(A, D)[:2]
    s_AB = s_AB / 2
    n = int(s_AB / 1000)
    ds = s_AB / n
    Phi_A, Lambda_A = [radians(i) for i in A]
    A_ab = radians(A_ab)
    for i in range(n):
        M = a*(1-e2)/(sqrt((1-e2*sin(Phi_A)**2)**3))
        N = a / (sqrt(1-e2*(sin(Phi_A)**2)))

        Phi_przyrost = ds * cos(A_ab) / M
        Az_przyrost = ds * sin(A_ab) * tan(Phi_A) / N

        mid_phi = Phi_A + 1/2*Phi_przyrost
        mid_az = A_ab + 1/2*Az_przyrost

        M = a*(1-e2)/(sqrt((1-e2*sin(mid_phi)**2)**3))
        N = a / (sqrt(1 - e2 * (sin(mid_phi) ** 2)))

        Phi_przyrost = ds*cos(mid_az)/M
        Lambda_przyrost = ds*sin(mid_az)/(N*cos(mid_phi))
        Az_przyrost = sin(mid_az)*tan(mid_phi)*ds/N

        Phi_A = Phi_A + Phi_przyrost
        Lambda_A = Lambda_A + Lambda_przyrost
        A_ab = A_ab + Az_przyrost

    return [degrees(Phi_A), degrees(Lambda_A)], degrees(A_ab)


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


def azymuty_pomiedzy_punktami():
    Azymut, odwrotny = vincent(kivioji()[0], Mean)[1:]
    return [Azymut - degrees(pi), odwrotny - degrees(pi)]


def surface_area(A, B):
    for i in range(0, 2):
        A[i] = radians(A[i])
        B[i] = radians(B[i])

    e = sqrt(e2)
    Phi_A = sin(A[0])/(1 - e2*(sin(A[0])**2)) + log((1+e*sin(A[0]))/(1-e*sin(A[0])))/(2*e)
    Phi_B = sin(B[0])/(1 - e2*(sin(B[0])**2)) + log((1+e*sin(B[0]))/(1-e*sin(B[0])))/(2*e)

    b2 = (a * sqrt(1 - e2))**2
    area = b2*(B[1] - A[1])/2*(Phi_A - Phi_B)
    return round(area, 6)

print(f"")
print(f"Współrzędne punktu średniej szerokości: ")
print(f"")
print(f"   phi={angles_formatting((Mean[0]))}, lambda={angles_formatting((Mean[1]))}")
print(f"")
print(f"Azymut pomiędzy punktami AD: ")
print(f"")
print(f"   {angles_formatting(vincent(A, D)[1])}")
print(f"")
print(f"Azymut odwrotny AD: ")
print(f"")
print(f"   {angles_formatting(vincent(A, D)[2])}")
print(f"")
print(f"Współrzędne punktu środkowego: ")
print(f"Phi:{angles_formatting(kivioji()[0][0])}    lambda:{angles_formatting(kivioji()[0][1])}    Azymut:{angles_formatting(kivioji()[1])}")
print(f"")
print(f"Odleglosc miedzy punktem średniej szerokości, a środkowym: ")
print(f"")
print(f"{round(vincent(Mean, kivioji()[0])[0], 3)}m")
print(f"")
print(f"Azumyt:  ")
print(f"")
print(f"    {angles_formatting(azymuty_pomiedzy_punktami()[0])}")
print(f"")
print(f"Azymut odwrotny:   ")
print(f"")
print(f"    {angles_formatting(azymuty_pomiedzy_punktami()[1])}")
print(f"")
print(f"Pole powierzchni czworokąta: {surface_area(A, D)}m^2")
