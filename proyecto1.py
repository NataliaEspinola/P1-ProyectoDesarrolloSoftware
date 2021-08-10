import sys
import io
import math

DIAMETRO = 1
RADIO = DIAMETRO / 2
INICIO = 0
CM = 0.5

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
        self.urm = None
        self.cd = None

    def calculate_ufz_urt_urm_cd(self, Taus):
        #ufz
        condition = 73 * math.sqrt(Taus)
        if condition < 5:
            self.ufz = 2.5 * math.log(condition * self.z) + 5.5
        elif condition >= 5 and condition < 70:
            self.ufz = (2.5 * math.log(condition * self.z) + 5.5) - (2.5 * math.log(1 + 0.3 * condition))
        else:
            self.ufz = 2.5 * math.log(30 * self.z)
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
        #por hacer
        pass
    def masaVirtual(self):
        #por hacer
        pass
    def lift(self):
        #por hacer
        pass

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
    for i in range(len(sys.argv)):
        if i == 1:
            with open(sys.argv[i]) as f:
                prm.T, prm.dt = map(float, f.readline().split())
                prm.theta, prm.R, prm.Taus, prm.CL = map(float, f.readline().split())
                for line in f:
                    x, y, z, u, v, w = map(float, line.split())
                    ptc = particula(x, y, z, u, v, w)
                    ptc.calculate_ufz_urt_urm_cd(prm.Taus)
                    particulas.append(ptc)
#////////////////FIN Lectura de txt/////////////////////////
    print(particulas)
    #INICIO DE CALCULOS
    try:
        print(particulas[0].cd)
    except Exception as e:
        print(f"{e}")