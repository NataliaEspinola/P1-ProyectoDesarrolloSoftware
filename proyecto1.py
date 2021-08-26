import sys
import io
from math import sqrt, pow, pi, log, pow, atan, tan, cos, sin
import random
import time
import multiprocessing as mp
from multiprocessing.pool import Pool

DIAMETRO = 1
RADIO = DIAMETRO / 2
INICIO = 0
CM = 0.5
radiandes_de_75_grados = 5 * pi / 12
radianes_de_10_grados = 0.174533
COUNTER = 0
OUT = ''

T = 0
dt = 0
theta = 0
R = 0
Taus = 0
CL = 0

'''
Fdrx, Fdry, Fdrz    FUERZAS DE ARRASTRE
Fswx, Fswz          PESO SUMERGIDO
Fvmx                MASA VIRTUAL
Flfz                FUERZA DE ELEVACION

'''

class particula:
    def __init__(self, x, y, z, u, v, w, Fdrx = None, Fdry = None, Fdrz = None, Fswx = None, Fswz = None, Fvmx = None, Flfz = None) -> None:
        self.x = x
        self.y = y
        self.z = z
        self.u = u
        self.v = v
        self.w = w
        self.urt = None
        self.ufz = None
        self.uftop = None
        self.ufbot = None
        self.urm = None
        self.cd = None
        self.ur2t = None
        self.ur2b = None
        self.saltos = 0
        self.max_z = z
        self.z_por_salto = []

def calculate_ufz_urt_urm_cd_uftop_ufbot_ur2t_ur2b(p, Taus):
    # ufz
    condition = 73 * sqrt(Taus)
    t = p.z - CM
    if t <= 0:
        t = 0.001
    if condition < 5:
        p.ufz = 2.5 * log(condition * p.z) + 5.5
        p.uftop = 2.5 * log(condition * (p.z + CM)) + 5.5
        p.ufbot = 2.5 * log(condition * t) + 5.5
    elif condition >= 5 and condition < 70:
        p.ufz = (2.5 * log(condition * p.z) + 5.5) - (2.5 * log(1 + 0.3 * condition))
        p.uftop = (2.5 * log(condition * (p.z + CM)) + 5.5) - (2.5 * log(1 + 0.3 * condition))
        p.ufbot = (2.5 * log(condition * t) + 5.5) - (2.5 * log(1 + 0.3 * condition))
    else:
        p.ufz = 2.5 * log(30 * p.z)
        p.uftop = 2.5 * log(30 * (p.z + CM))
        p.ufbot = 2.5 * log(30 * t)
    # urt
    p.urt = p.u - p.ufz
    # urm
    p.urm = sqrt(pow(p.urt, 2) + pow(p.v, 2) + pow(p.w, 2))
    # cd
    rep = p.urm * condition
    p.cd = 24 / (rep * (1 + 0.15 * sqrt(rep) + 0.017 * rep) - (
                0.208 / (1 + pow(10, 4) * pow(rep, -0.5))))
    # ur2t
    p.ur2t = pow((p.u - p.uftop), 2) + pow(p.v, 2) + pow(p.w, 2)
    # ur2b
    p.ur2b = pow((p.u - p.ufbot), 2) + pow(p.v, 2) + pow(p.w, 2)
    return p

def drag(p, R):
    comun = -0.75 * (1 / (1 + R + CM)) * p.cd * p.urm
    p.Fdrx = comun * p.urt
    p.Fdry = comun * p.v
    p.Fdrz = comun * p.w
    return p

def pesoSumergido(p, theta, Taus, R):
    # Comun
    comun = (1 / (1 + R + CM)) * (1 / Taus)
    # Fswx
    p.Fswx = sin(theta) * comun
    # Fswz
    p.Fswz = cos(theta) * -comun
    return p

def masaVirtual(p, R):
    # Fvmx
    p.Fvmx = (CM / (1 + R + CM)) * p.w * (2.5 / p.z)
    return p

def lift(p, R, CL):
    # Flfz
    p.Flfz = 0.75 * (1 / (1 + R + CM)) * CL * (p.ur2t - p.ur2b)
    return p

def efecto_choque(p):  # da nuevas velocidades
    # w luego del rebote
    new_w = -p.w
    p.w = new_w
    # u luego del rebote
    e = random.uniform(0.0, radianes_de_10_grados)  # Random float:  0.0 <= x <= 10.0
    alpha = atan(new_w / p.u)
    while alpha >= radiandes_de_75_grados:  # compara en radianes
        e = random.uniform(0.0, radianes_de_10_grados)
        alpha = atan(new_w / p.u)
    new_u = new_w / tan(alpha + e)
    p.u = new_u
    # v luego del rebote
    angulo_para_Y = random.uniform(-radianes_de_10_grados,
                                    radianes_de_10_grados)  # Random float:  -10.0 <= x <= 10.0
    new_v = new_u * tan(angulo_para_Y)
    p.v = new_v
    return p

def new_u_v_w(p, dt):
    # u
    p.u = p.u + (dt * (p.Fdrx + p.Fswx + p.Fvmx))
    # v
    p.v = p.v + (dt * (p.Fdry))
    # w
    p.w = p.w + (dt * (p.Fdrz + p.Flfz + p.Fswz))
    return p

def new_x_y_z(p, dt):
    # x
    p.x = p.x + p.u * dt
    # y
    p.y = p.y + p.v * dt
    # z
    p.z = p.z + p.w * dt
    return p
def new_pos_choque(p, dt):
    # x
    p.x = p.x + p.u * dt
    # y
    p.y = p.y + p.v * dt
    # z
    p.z = 0.501
    return p


if __name__ == "__main__":
    inicio = time.time()
    particulas = []
    largo = 0
#////////////////INICIO Lectura de txt/////////////////////////    
    for i in range(len(sys.argv)):
        if i == 1:
            with open(sys.argv[i]) as f:
                OUT = sys.argv[i].split('.', 1)[0]
                T, dt = map(float, f.readline().split())
                theta, R, Taus, CL = map(float, f.readline().split())
                for p in f:
                    x, y, z, u, v, w = map(float, p.split())
                    ptc = particula(x, y, z, u, v, w)
                    ptc = calculate_ufz_urt_urm_cd_uftop_ufbot_ur2t_ur2b(ptc,Taus)
                    ptc = pesoSumergido(ptc,theta, Taus, R)
                    ptc = masaVirtual(ptc,R)
                    ptc = drag(ptc,R)
                    ptc = lift(ptc,R, CL)
                    largo += 1
                    particulas.append(ptc)
#////////////////FIN Lectura de txt/////////////////////////
    #INICIO DE CALCULOS
    while COUNTER < T:
        COUNTER += dt
        with mp.Pool(processes=4) as pool:                                          
            for t in range(largo):                                                   
                particulas = pool.starmap(calculate_ufz_urt_urm_cd_uftop_ufbot_ur2t_ur2b, particulas)
        for i in range(len(particulas)):
            # se guarda la altura anterior
            last_z = particulas[i].z
            # ver rebote y asignar un +1 al salto si paso
            if particulas[i].z < 0.501:
                particulas[i] = efecto_choque(particulas[i])
                particulas[i].saltos += 1
                particulas[i] = new_pos_choque(particulas[i],dt)
            else:
                # Nueva vel
                particulas[i] = new_u_v_w(particulas[i],dt)
                # pos
                particulas[i] = new_x_y_z(particulas[i],dt)
                # comparar altura anterior con actual para guardar altura por salto
                if particulas[i].saltos > 0:
                    if last_z > particulas[i].z:
                        particulas[i].z_por_salto.append(last_z)
            # fuerzas
            particulas[i] = calculate_ufz_urt_urm_cd_uftop_ufbot_ur2t_ur2b(particulas[i],Taus)
            particulas[i] = masaVirtual(particulas[i],R)
            particulas[i] = drag(particulas[i],R)
            particulas[i] = lift(particulas[i],R, CL)
            # ver z max
            if particulas[i].z > particulas[i].max_z:
                particulas[i].max_z = particulas[i].z


    salida = OUT + ".out"
    with open(salida, 'w') as f:
        a = 0
        for p in particulas:
            if len(p.z_por_salto) > 0:
                f.write("".join([str(p.x), " ", str(p.y), " ", str(p.z), " ", str(p.saltos), " ", str(p.max_z), " ",
                                 str(sum(p.z_por_salto) / len(p.z_por_salto)), "\n"]))
            else:
                f.write("".join([str(p.x), " ", str(p.y), " ", str(p.z), " ", str(p.saltos), " ", str(p.max_z), " 0",
                                 "\n"]))
            # print(
            #     f"Particula {a}\n        x: {p.x}, y: {p.y}, z: {p.z}\n     Fdrx: {p.Fdrx}, Fswx: {p.Fswx}, Fvmx: {p.Fvmx}\n     Fdry: {p.Fdry}\n     Fdrz: {p.Fdrz}, Fwsz: {p.Fswz}, Flfz: {p.Flfz}")
            # print("ut: ", p.u, ", vt: ", p.v, ",wt: ", p.w)
            # a += 1
    fin = time.time()
    print("Tiempo: ", fin - inicio)