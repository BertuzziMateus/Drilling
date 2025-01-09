
import numpy as np

class DataSet:

    def __init__(self,
                P0 : tuple[float,float], # X and Y coordenates from the initial point
                P3 : tuple[float,float], # X and Y coordenates objective
                ro_fluid: float, #Ro of the fluid
                ro_command: float, #Ro of the command pipe
                ro_drillpipe: float, #Ro of the drill pipe
                ro_heavypipe: float, #Ro of the heavy pipe
                diameters_command: tuple[float,float], # Informations about the command pipe [De,Di]
                diameters_drillpipe: tuple[float,float], # Informations about the drill pipe [De,Di]
                diameters_heavypipe: tuple[float,float], # Informations about the heavy pipe [De,Di]
                lp : float, # len of heavypipe
                µ: float, #Coefficient of friction
                z: float, # Rock reaction 
                max: float,
                radius: tuple[float,float],
                ) -> None:
        area_command = ( np.pi/4)*(diameters_command[0]**2  - diameters_command[1]**2)
        area_drill =(np.pi/4)*(diameters_drillpipe[0]**2  - diameters_drillpipe[1]**2)
        area_heavy =(np.pi/4)*(diameters_heavypipe[0]**2  - diameters_heavypipe[1]**2)
        self.P0 = P0 
        self.P3 = P3
        self.g = 9.81
        self.ro_fluid = ro_fluid #kg / m^3 
        self.ro_command = ro_command #kg / m^3 
        self.ro_drillpipe = ro_drillpipe #kg / m^3 
        self.ro_heavypipe = ro_heavypipe #kg / m^3 
        self.µ = µ
        self.z = z
        self.area_command = area_command
        self.area_drill = area_drill
        self.area_heavy = area_heavy
        self.lambd_command = 148.816 #Kg/m
        self.lambd_drill = 24.55471 #Kg/m
        self.lambd_heavy = 74.4082 #Kg/m
        self.d_ext_drill = diameters_drillpipe[0]
        self.d_int_drill = diameters_drillpipe[1]
        self.d_ext_command = diameters_command[0]
        self.d_int_command = diameters_command[1]
        self.d_ext_heavy = diameters_heavypipe[0]
        self.d_int_heavy = diameters_heavypipe[1]
        self.lp = lp
        self.max = max # m
        self.min_radius = radius[0]
        self.max_radius = radius[1]

        pass