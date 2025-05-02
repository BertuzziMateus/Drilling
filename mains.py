import Auxiliaries as ax
from Datas import *
import numpy as np

def minimal_tension(Data) -> list:
    t1_data = []
    l1_data = []
    R_data = []

    l1 = 100
    condition =  False 
    while condition == False:

        t1_data_R = []
        R_valor = []

        for R in np.arange(Data.min_radius,Data.max_radius+50,50):
            if ax.theta(Data,l1,R)*(180/np.pi) > 52:
                pass
            else: 
                up_t1,*_ = ax.up_tension(Data,l1,R)
                t1_data_R.append(up_t1)
                R_valor.append(R)
        if len(t1_data_R) != 0 :
            key_r = np.where(t1_data_R == min(t1_data_R))[0][0]
            t1_data.append(t1_data_R[key_r])
            l1_data.append(l1)
            R_data.append(R_valor[key_r])

        if round(l1,1) == Data.max:
            condition = True

        l1 += 10

    key = np.where(t1_data == min(t1_data))[0][0]
    l1_value = l1_data[key]
    R_value = R_data[key]




    return [l1_value,R_value]

def minimal_torque(Data) -> list:
    torque_data = []
    l1_data = []
    R_data = []

    l1 = 100

    condition =  False 

    while condition == False:

        torque_data_R = []
        R_valor = []

        for R in np.arange(Data.min_radius,Data.max_radius+50,50):
            if ax.theta(Data,l1,R)*(180/np.pi) > 52:
                pass
            else: 
                t1,t2,t3,torque = ax.down_tension(Data,l1,R)
                torque_data_R.append(torque)
                R_valor.append(R)

        if len(torque_data_R) != 0 :
            key_r = np.where(torque_data_R == min(torque_data_R))[0][0]
            torque_data.append(torque_data_R[key_r])
            l1_data.append(l1)
            R_data.append(R_valor[key_r])

        if round(l1,1) == Data.max:
            condition = True

        l1 += 10

    key = np.where(torque_data == min(torque_data))[0][0]
    l1_value = l1_data[key]
    R_value = R_data[key]

    return [l1_value,R_value]


def time_drilling(Data,l1,R,Lithology) -> float:

    l1,l2,l3 = ax.lenght(Data,l1,R)

    total_length = l1 + l2 + l3

    lithology_mesh = Lithology_mesh(Data,Lithology)
    mesh = lithology_mesh.mesh
    mesh_lithology = lithology_mesh.mesh_lithology

    time_in_hours = 0.0

    for i in range(len(mesh)):
        lithology  = mesh_lithology[i]
        height = mesh[i]
        rop = Rop_rocks[lithology]
        time = height / rop
        time_in_hours += time

    print(f'Time in hours: {time_in_hours:.2f}')
    print(f'Time in days: {(time_in_hours/24):.2f}')

    return time_in_hours


def time_drilling_real(Data,l1,R,Lithology) -> float:

    l1,l2,l3 = ax.lenght(Data,l1,R)

    total_length = l1 + l2 + l3

    lithology_mesh = Lithology_mesh(Data,Lithology)
    mesh = lithology_mesh.mesh
    mesh_lithology = lithology_mesh.mesh_lithology

    dx = 0.0

    while round(dx)!= round(total_length):

        


        dx += 1



