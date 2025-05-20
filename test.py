from Mains import *
from Results import *


Data = Well_data(P0=(0, 0),  # X , Y
                 P3=(1000, 3000),  # X , Y
                 ro_fluid=1737.5,  # kg / m^3
                 ro_command=8000,  # kg / m^3
                 ro_drillpipe=8000,  # kg / m^3
                 ro_heavypipe=8000,  # kg / m^3
                 diameters_command=(0.2032, 0.0714502),  # m 8 e 2.813
                 diameters_drillpipe=(0.127,  0.1086104),  # m 5 e 4.276
                 diameters_heavypipe=(0.1524,  0.1016),  # m 6 e 4
                 µ=0.23,
                 z=(5000*8) * 4.44822,  # Newton
                 lp=36,  # m
                 max=2300,
                 radius=(100, 600)
                 )

litologia = [
    ['Sandstone', 500],
    ['Limestone', 400],
    ['Evaporite', 1000],
    ['Shale', 1000],
    ['Sandstone', 100],
]

results(Data,litologia)
Results_for_drilling(Data)



