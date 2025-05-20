
import numpy as np

class Well_data:

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
        self.lambd_command = ro_command * area_command
        self.lambd_drill = ro_drillpipe * area_drill
        self.lambd_heavy = ro_heavypipe * area_heavy
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



class Lithology_mesh:
    def __init__(self, well_data, list_lithology) -> None:
        self.list_lithology = list_lithology
        
        def make_mesh(litho_list):
            mesh_height = []
            mesh_lithology = []
            for element in litho_list:
                mesh_height.append(element[1])
                mesh_lithology.append(element[0])
            return mesh_height, mesh_lithology
        
        
        self.mesh_height, self.mesh_lithology = make_mesh(list_lithology)

        mesh = []
        for element in self.mesh_height:
            mesh.append(element)
        
        self.mesh = mesh


Colors_rocks = {
    'Sandstone': '#ffff40',
    'Limestone': '#00efef',
    #'Dolomite': '#00bfcf',
    'Evaporite': '#799fbf',
    'Shale': '#40ff00',
}





class Sandstone:
    def __init__(self,
                 ) -> None:
        a = 3.14e-7
        b = 94.58
        rop = a*b
        self.rop = 10 #m/h
    pass
pass

class Limestone:
    def __init__(self,
                 ) -> None:
        a = 2.33e-7
        b = 92.31
        rop = a*b
        self.rop = 10
    pass
pass

class Evaporite:
    def __init__(self,
                 ) -> None:
        
        a =1.74e-7
        b = 93.85
        rop = a*b
        self.rop = 10
    pass
pass
class Shale:
    def __init__(self,
                 ) -> None:
        a = 3.54e-7
        b = 76.14
        rop = a*b
        self.rop = 10
    pass
pass


Rop_rocks = {
    'Sandstone': Sandstone().rop, #ft/h
    'Limestone': Limestone().rop, #ft/h
    #'Dolomite': 30, #ft/h
    'Evaporite': Evaporite().rop, #ft/h
    'Shale': Shale().rop, #ft/h
}


        

    

    



        









# class Rock:

#     class Sandstone:

#         def __init__(self,
#                 color : str,
#                 ):
#             self.rop = 20
#             self.color = color
#             return None
#         pass

#     class Limestone:   
#         def Limestone(self,
#             color: str = '48d1cc',
#             ):
#             self.rop  = 30
#             return None

#     class Dolomite:  
#         def Dolomite(self,
#             color: str = '53c5ec',
#             ):
#             self.rop = 15
#             return None
        
#     class Evaporite:
#         def Evaporite(self,
#             color: str = '799fbf',
#             ):
#             self.rop = 10
#             return None

#         pass
    
