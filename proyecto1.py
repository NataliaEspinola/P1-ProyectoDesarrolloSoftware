import sys
import io
import math
import random

DIAMETRO = 1
RADIO = DIAMETRO / 2
INICIO = 0
CM = 0.5
radiandes_de_75_grados = 5*math.pi/12
radianes_de_10_grados = 0.174533
COUNTER = 0
OUT = ''

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
        self.promedio_alturas = None
        self.z_por_salto = []

    def calculate_ufz_urt_urm_cd_uftop_ufbot(self, Taus):
        #ufz
        condition = 73 * math.sqrt(Taus)
        t = self.z - CM
        if t <= 0:
            t = 0.01
        if condition < 5:
            self.ufz = 2.5 * math.log(condition * self.z) + 5.5
            self.uftop = 2.5 * math.log(condition * (t)) + 5.5
            self.ufbot = 2.5 * math.log(condition * t) + 5.5
        elif condition >= 5 and condition < 70:
            self.ufz = (2.5 * math.log(condition * self.z) + 5.5) - (2.5 * math.log(1 + 0.3 * condition))
            self.uftop = (2.5 * math.log(condition * (self.z + CM)) + 5.5) - (2.5 * math.log(1 + 0.3 * condition))
            self.ufbot = (2.5 * math.log(condition * t) + 5.5) - (2.5 * math.log(1 + 0.3 * condition))
        else:
            self.ufz = 2.5 * math.log(30 * self.z)
            self.uftop = 2.5 * math.log(30 * (self.z + CM))
            self.ufbot = 2.5 * math.log(30 * t)
        #urt
        self.urt = self.u - self.ufz
        #urm
        self.urm = math.sqrt(math.pow(self.urt,2) + math.pow(self.v,2) + math.pow(self.w,2))
        #cd
        rep = self.urm * condition
        self.cd = 24 / ( rep * (1 + 0.15 * math.sqrt(rep) + 0.017 * rep) - (0.208/(1 + math.pow(10,4) * math.pow(rep,-0.5))))

    def drag(self, R):
        comun = -0.75 * (1 / 1+R+CM) * self.cd * self.urm
        self.Fdrx = comun * self.urt
        self.Fdry = comun * self.v
        self.Fdrz = comun * self.w

    def pesoSumergido(self, theta, Taus, R):
        #Comun
        comun = ( 1 / 1 + R + CM ) * ( 1 / Taus )
        #Fswx
        self.Fswx = math.sin(theta) * comun
        #Fswz
        self.Fswz = math.cos(theta) * -comun

    def masaVirtual(self, R):
        #Fvmx
        self.Fvmx = ( CM / 1 + R + CM ) * self.w * ( 2.5 / self.z)

    def lift(self, R, CL):
        #Flfz
        self.Flfz = 0.75 * (1 / 1 + R + CM) * CL * (self.ur2t + self.ur2b)

    def calculate_ur2t_ur2b(self):
        #ur2t
        self.ur2t = math.pow(self.u - self.uftop, 2) + math.pow(self.v, 2) + math.pow(self.w, 2)
        #ur2b
        self.ur2b = math.pow(self.u - self.ufbot, 2) + math.pow(self.v, 2) + math.pow(self.w, 2)

    def efecto_choque(self): # da nuevas velocidades
        # w luego del rebote
        new_w = -self.w
        self.w = new_w #revisar
        # u luego del rebote
        e = random.uniform(0.0, radianes_de_10_grados) # Random float:  0.0 <= x <= 10.0
        alpha = math.atan(new_w/self.u)
        while alpha >= radiandes_de_75_grados: # compara en radianes
            e = random.uniform(0.0, radianes_de_10_grados)
            alpha = math.atan(new_w / self.u)
        new_u = new_w/math.tan(alpha+e)
        self.u = new_u #revisar
        # v luego del rebote
        angulo_para_Y = random.uniform(-radianes_de_10_grados, radianes_de_10_grados) # Random float:  -10.0 <= x <= 10.0
        new_v = new_u * math.tan(angulo_para_Y)
        self.v = new_v #revisar
    
    def new_u_v_w(self,dt):
        #u
        self.u = self.u + ( dt * (self.Fdrx + self.Fswx + self.Fvmx))
        #v
        self.v = self.v + ( dt * (self.Fdry))
        #w
        self.w = self.w + ( dt * (self.Fdrz + self.Flfz + self.Fswz))
    def new_x_y_z(self, dt):
        #x
        self.x = self.x + self.u * dt
        #y
        self.y = self.y + self.v * dt
        #z
        self.z = self.z + self.w * dt


class parametros:
    def __init__(self, T = None, dt = None, theta = None, R = None, Taus = None, CL = None) -> None:
        self.T = T
        self.dt = dt
        self.theta = theta
        self.R = R
        self.Taus = Taus
        self.CL = CL

if __name__ == "__main__":
    particulas = []
    prm = parametros()
#////////////////INICIO Lectura de txt/////////////////////////
#/////////////Debe recibir el txt como argumento en consola/////
    for i in range(len(sys.argv)):
        if i == 1:
            try:
                with open(sys.argv[i]) as f:
                    OUT = sys.argv[i].strip('.txt')
                    prm.T, prm.dt = map(float, f.readline().split())
                    prm.theta, prm.R, prm.Taus, prm.CL = map(float, f.readline().split())
                    for line in f:
                        x, y, z, u, v, w = map(float, line.split())
                        ptc = particula(x, y, z, u, v, w)
                        ptc.calculate_ufz_urt_urm_cd_uftop_ufbot(prm.Taus)
                        ptc.pesoSumergido(prm.theta, prm.Taus, prm.R)
                        ptc.masaVirtual(prm.R)
                        ptc.drag(prm.R)
                        ptc.calculate_ur2t_ur2b()
                        ptc.lift(prm.R, prm.CL)
                        particulas.append(ptc)
            except Exception as e:
                print('{e}')
#////////////////FIN Lectura de txt/////////////////////////
    #INICIO DE CALCULOS
    try:
        #SIMULACION
        while COUNTER < prm.T:
            COUNTER += prm.dt
            for i in range(len(particulas)):
                # ver rebote y asignar un +1 al salto si paso
                if particulas[i].z < 0.501:
                    particulas[i].efecto_choque()
                    particulas[i].saltos += 1
                    particulas[i].new_x_y_z(prm.dt)
                    particulas[i].z_por_salto.append(particulas[i].max_z)
                else:
                    #Nueva vel
                    particulas[i].new_u_v_w(prm.dt)
                    #pos
                    particulas[i].new_x_y_z(prm.dt)
                #fuerzas
                particulas[i].calculate_ufz_urt_urm_cd_uftop_ufbot(prm.Taus)
                particulas[i].masaVirtual(prm.R)
                particulas[i].drag(prm.R)
                particulas[i].calculate_ur2t_ur2b()
                particulas[i].lift(prm.R, prm.CL)
                #ver z max
                if particulas[i].z > particulas[i].max_z:
                    particulas[i].max_z = particulas[i].z

        salida = OUT+".out"
        with open(salida,'w') as f:
            for line in particulas:
                f.write(str(line.x)+' '+str(line.y)+' '+str(line.z)+' '+str(line.saltos)+' '+str(line.max_z)+' '+str(sum(line.z_por_salto)/len(line.z_por_salto))+'\n')
    except Exception as e:
        print(f"{e}")