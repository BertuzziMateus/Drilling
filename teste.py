from Drilling import *

import time


Data = DataSet(P0 = (0, 0), # X , Y
                P3 = (1000, 3000),  # X , Y
                ro_fluid = 1737.5 , #kg / m^3 
                ro_command = 6712.96, #kg / m^3 
                ro_drillpipe = 10201.97, #kg / m^3 
                ro_heavypipe = 9323.58, #kg / m^3 
                diameters_command = (0.2032, 0.1143) , # m 8 e 4.5
                diameters_drillpipe = (0.127, 0.1143), # m 5 e 4.5
                diameters_heavypipe = (0.1524, 0.1143), # m 6 e 4.5
                µ = 0.23, # 
                z = (5000*8) * 4.44822 , # Newton
                lp = 36, # m 
                max = 2300,
                radius= (100,600)
                )

l1,r = minimal_tension(Data)
# l1_1,r_1 = minimal_torque(Data)
drilling_informations_table(Data)
drilling_draw(Data)
#tension_in_radius(Data,l1)
tension_in_section1(Data,r)
tension_graphic(Data)
