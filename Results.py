import Auxiliaries as ax
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Mains import *


def results(Data,Lithology):

    l1_tension,R1_tension = minimal_tension(Data)
    l1_torque,R1_torque = minimal_torque(Data)

    def drilling_draw(Data,Lithology) -> None :

        x1, y1 = ax.points_coordinates(Data, l1_tension,R1_tension)
        x2, y2 = ax.points_coordinates(Data, l1_torque,R1_torque)

        lithology_mesh = Lithology_mesh(Data,Lithology)
        mesh = lithology_mesh.mesh
        mesh_lithology = lithology_mesh.mesh_lithology

        colors = []


        for element in mesh_lithology:
            colors.append(element)

        custom_style = {
            'font.size': 16,  # Tamanho adequado para leitura de gráficos
            'axes.labelsize': 18,  # Tamanho dos rótulos dos eixos
            'axes.titlesize': 14,  # Tamanho do título do gráfico
            'axes.linewidth': 1.5,  # Espessura das bordas dos gráficos
            'xtick.labelsize': 18,  # Tamanho do texto dos ticks no eixo x
            'ytick.labelsize': 18,  # Tamanho do texto dos ticks no eixo y
            'lines.linewidth': 2,  # Espessura das linhas dos gráficos
            'lines.markersize': 6,  # Tamanho dos marcadores
            'legend.fontsize': 18,  # Tamanho da legenda
            'legend.frameon': False,  # Remove a moldura ao redor da legenda
            'legend.loc': 'best',  # Melhor posição automática para a legenda
            'figure.figsize': (4, 3),  # Tamanho padrão da figura (polegadas)
            'savefig.dpi': 1000,  # Alta resolução para exportação (publicação)
            'savefig.bbox': 'tight',  # Salva a imagem sem cortar parte do gráfico
            }

        plt.rcParams.update(custom_style)
        plt.axes().set_aspect('equal')
        plt.plot(x1, y1, color="k") 
        plt.gca().invert_yaxis()
        plt.title("Type 1 well trajectory for the minimal tension")
        plt.xlabel("Distance ($m$)")
        plt.ylabel("Depth ($m$)")
        plt.grid(alpha=0.7)
        mesh_old = 0
        sum_mesh = mesh[0]
        for i in range(len(mesh)):
            plt.axhspan(mesh_old, sum_mesh, color=Colors_rocks[mesh_lithology[i]], alpha=0.5, label=mesh_lithology[i])
            mesh_old += mesh[i]
            if i == len(mesh)-1:
                sum_mesh += mesh[i]
            else:
                sum_mesh += mesh[i+1]
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
        plt.show()

        plt.axes().set_aspect('equal')
        plt.plot(x2, y2, color="red", linestyle="--")
        plt.gca().invert_yaxis()
        plt.title("Type 1 well trajectory for the minimal Torque")
        plt.xlabel("Distance ($m$)")
        plt.ylabel("Depth ($m$)")
        plt.grid(alpha=0.7)
        mesh_old = 0
        sum_mesh = mesh[0]
        for i in range(len(mesh)):
            plt.axhspan(mesh_old, sum_mesh, color=Colors_rocks[mesh_lithology[i]], alpha=0.5, label=mesh_lithology[i])
            mesh_old += mesh[i]
            if i == len(mesh)-1:
                sum_mesh += mesh[i]
            else:
                sum_mesh += mesh[i+1]
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
        plt.show()


        pass

    def drilling_informations(Data)->list:

        up_tensions = ax.up_tension(Data,l1_tension,R1_tension)
        down_tensions = ax.down_tension(Data,l1_tension,R1_tension)
        angle = ax.theta(Data,l1_tension,R1_tension)*180/np.pi
        neutral_line = ax.Nl(Data,l1_tension,R1_tension)
        lenghts = ax.lenght(Data,l1_tension,R1_tension)
        lenght_command = ax.buckling(Data,l1_tension,R1_tension)

        tension_results = [up_tensions,down_tensions,angle,neutral_line,lenghts,lenght_command,l1_tension,R1_tension]

        up_tensions = ax.up_tension(Data,l1_torque,R1_torque)
        down_tensions = ax.down_tension(Data,l1_torque,R1_torque)
        angle = ax.theta(Data,l1_torque,R1_torque)*180/np.pi
        neutral_line = ax.Nl(Data,l1_torque,R1_torque)
        lenghts = ax.lenght(Data,l1_torque,R1_torque)
        lenght_command = ax.buckling(Data,l1_torque,R1_torque)

        torque_results = [up_tensions,down_tensions,angle,neutral_line,lenghts,lenght_command,l1_torque,R1_torque]

        return [tension_results,torque_results]
    
    def drilling_informations_table(data):
        results = drilling_informations(data)
        for i in range(len(results)):
            up_tensions,down_tensions,angle,neutral_line,lenghts,lenght_command,l1,R = results[i]
            l1,l2,l3 = lenghts
            t1_up ,t2_up,t3_up = up_tensions
            t1_down,t2_down,t3_down,torque = down_tensions
            data = [l1,l2,l3,R,t1_up,t2_up,t3_up,t1_down,t2_down,t3_down,torque,angle,neutral_line,lenght_command ]
            data = np.round(data,2)
            index_labels = [
                "L1 (m):", "L2 (m):", "L3 (m):","radius (m):",
                "Up tension L1 (N/m²):", "Up tension L2 (N/m²):", "Up tension L3 (N/m²):",
                "Down tension L1 (N/m²):", "Down tension L2 (N/m²):", "Down tension L3 (N/m²):",
                "Torque (Nm):", "Angle (°):", "Neutral line (m):", "Length command (m):"]
            
            if i == 0:
                print(f"\n--- Result table for minimal tension ---")
            else:
                print(f"\n--- Result table for minimal torque ---")
            a = pd.DataFrame(data,columns=[''],index=index_labels)
            print(a)
            print('')


        pass


    

    def time_drilling_real(Data, l1, R, Lithology) -> float:

        l1, l2, l3 = ax.lenght(Data, l1, R)
        angle = ax.theta(Data, l1, R)
        lithology_mesh = Lithology_mesh(Data, Lithology)

        mesh = lithology_mesh.mesh
        mesh_lithology = lithology_mesh.mesh_lithology
        grid_dict = dict()
        sum_l = 0
        indictor = 0

        for element in mesh:
            for _ in range(element):
                grid_dict[f'{sum_l}'] = f'{mesh_lithology[indictor]}'
                sum_l += 1
            grid_dict[f'{sum_l}'] = f'{mesh_lithology[indictor]}'
            indictor += 1

        l3_v = (l3*np.cos(angle))
        l2_v = Data.P3[1] - l1 - l3_v


        l1_v_for_dict = np.arange(0, l1, 1)
        l2_v_for_dict = np.round(np.arange(l1, l1 + l2_v, 1))
        l3_v_for_dict = np.round(np.arange(l1+l2_v, l1+l2_v+l3_v, 1))

        time_l1 = 0
        for element in l1_v_for_dict:
            lithology_type = grid_dict[f'{element}']
            time_l1 += 1*Rop_rocks[lithology_type]

        l2_vec = np.round(np.linspace(l1, l1+l2, len(l2_v_for_dict)))

        count = 0
        time_l2 = 0
        for element in l2_v_for_dict:
            lithology_type = grid_dict[f'{int(element)}']

            if count == 0:
                part1 = l2_vec[count]
                part2 = l1
            else:
                part1 = l2_vec[count]
                part2 = l2_vec[count-1]

            x_l2 = part1 - part2
            time_l2 += x_l2*Rop_rocks[lithology_type]
            count += 1

        l3_vec = np.round(np.linspace(l1+l2, l1+l2+l3, len(l3_v_for_dict)))

        count = 0
        time_l3 = 0
        for element in l3_v_for_dict:

            if count == 0:
                part1 = l3_vec[count]
                part2 = l1+l2
            else:
                part1 = l3_vec[count]
                part2 = l3_vec[count-1]

            lithology_type = grid_dict[f'{int(element)}']
            x_l3 = part1 - part2
            time_l3 += x_l3*Rop_rocks[lithology_type]
            count += 1
            
        time_h = time_l1+time_l2+time_l3

        print(f'Time in hours: {time_h:.2f}')
        print(f'Time in days: {(time_h/24):.2f}')
        return time_h

    drilling_draw(Data,Lithology)
    drilling_informations_table(Data)
    print(f'Time Perfuration For the Minimal Tension:')
    time_drilling_real(Data,l1_tension,R1_tension,Lithology)
    print('\n')
    print(f'Time Perfuration For the Minimal Torque:')
    time_drilling_real(Data,l1_torque,R1_torque,Lithology)

    

    pass


def Results_for_drilling(Data):
        t1_data = []
        l1_data = []
        R_data = []
        torque_data = []

        l1 = 100

        condition =  False 

        while condition == False:

            t1_data_R = []
            R_valor = []
            torque_data_R = []

            for R in np.arange(Data.min_radius,Data.max_radius+50,50):
                if ax.theta(Data,l1,R)*(180/np.pi) > 52:
                    pass
                else: 
                    up_t1,*_ = ax.up_tension(Data,l1,R)
                    t1,t2,t3,torque = ax.down_tension(Data,l1,R)
                    t1_data_R.append(up_t1)
                    torque_data_R.append(torque)
                    R_valor.append(R)

            if len(t1_data_R) != 0 :
                key_r = np.where(t1_data_R == min(t1_data_R))[0][0]
                key_r1 = np.where(torque_data_R == min(torque_data_R))[0][0]
                t1_data.append(t1_data_R[key_r])
                l1_data.append(l1)
                R_data.append(R_valor[key_r])
                torque_data.append(torque_data_R[key_r1])

            if round(l1,1) == Data.max:
                condition = True

            l1 += 10

        custom_style = {
                'font.size': 16,  # Tamanho adequado para leitura de gráficos
                'axes.labelsize': 18,  # Tamanho dos rótulos dos eixos
                'axes.titlesize': 14,  # Tamanho do título do gráfico
                'axes.linewidth': 1.5,  # Espessura das bordas dos gráficos
                'xtick.labelsize': 18,  # Tamanho do texto dos ticks no eixo x
                'ytick.labelsize': 18,  # Tamanho do texto dos ticks no eixo y
                'lines.linewidth': 2,  # Espessura das linhas dos gráficos
                'lines.markersize': 6,  # Tamanho dos marcadores
                'legend.fontsize': 18,  # Tamanho da legenda
                'legend.frameon': False,  # Remove a moldura ao redor da legenda
                'legend.loc': 'best',  # Melhor posição automática para a legenda
                'figure.figsize': (4, 3),  # Tamanho padrão da figura (polegadas)
                'savefig.dpi': 1000,  # Alta resolução para exportação (publicação)
                'savefig.bbox': 'tight',  # Salva a imagem sem cortar parte do gráfico
                }


        plt.rcParams.update(custom_style)

        plt.figure(figsize=(12, 6)) 
        plt.plot(l1_data,t1_data, label="Tension", color="blue",)
        plt.title("Tension behavior depending on the best configurations for section 1", fontsize=10)
        plt.xlabel("Length ($m$)")
        plt.ylabel("Tension ($N/m²$)")
        plt.legend()
        plt.margins(x=0.1, y=0.1)
        plt.grid(alpha=1,linewidth=2)
        plt.show()


        plt.plot(l1_data,torque_data, label="Torque", color="red", linestyle="--")
        plt.title("Torque behavior depending on the best configurations for section 1", fontsize=10)
        plt.xlabel("Length ($m$)")
        plt.ylabel("Torque ($N*m$)")
        plt.legend()
        plt.margins(x=0.1, y=0.1)
        plt.grid(alpha=1,linewidth=2)
        plt.show()

        pass











    
'''
    #Here we take the information obtained from the drilling informations function and using this to create a tension graph through radius.
    def tension_in_radius(Data,l1):

        torque_data = [] 
        t1_data_R = []
        R_valor = []

        for R in np.arange(Data.min_radius,Data.max_radius+50,50):
            if ax.theta(Data,l1,R)*(180/np.pi) > 52:
                pass
            else: 
                up_t1,*_ = ax.up_tension(Data,l1,R)
                t1,t2,t3,torque = ax.down_tension(Data,l1,R)
                t1_data_R.append(up_t1)
                torque_data.append(torque)
                R_valor.append(R)


        custom_style = {
            'font.size': 16,  # Tamanho adequado para leitura de gráficos
            'axes.labelsize': 18,  # Tamanho dos rótulos dos eixos
            'axes.titlesize': 14,  # Tamanho do título do gráfico
            'axes.linewidth': 1.5,  # Espessura das bordas dos gráficos
            'xtick.labelsize': 18,  # Tamanho do texto dos ticks no eixo x
            'ytick.labelsize': 18,  # Tamanho do texto dos ticks no eixo y
            'lines.linewidth': 2,  # Espessura das linhas dos gráficos
            'lines.markersize': 6,  # Tamanho dos marcadores
            'legend.fontsize': 18,  # Tamanho da legenda
            'legend.frameon': False,  # Remove a moldura ao redor da legenda
            'legend.loc': 'best',  # Melhor posição automática para a legenda
            'figure.figsize': (4, 3),  # Tamanho padrão da figura (polegadas)
            'savefig.dpi': 1000,  # Alta resolução para exportação (publicação)
            'savefig.bbox': 'tight',  # Salva a imagem sem cortar parte do gráfico
            }

        plt.rcParams.update(custom_style)

        plt.plot(R_valor,t1_data_R, label=f'Tension in section of ${l1_tension}m$', color="blue")
        plt.title("Tension behavior across the radius")
        plt.xlabel("Radius ($m$)")
        plt.ylabel("Tension ($N/m²$)")
        plt.legend()
        plt.grid(alpha=1,linewidth=2)
        plt.margins(x=0.1, y=0.1)
        plt.show()


        plt.plot(R_valor,torque_data, label=f'Torque in section of ${l1_tension}m$', color="red", linestyle="--")
        plt.title("Torque behavior across the radius")
        plt.xlabel("Radius ($m$)")
        plt.ylabel("Torque ($N*m$)")
        plt.legend()
        plt.grid(alpha=1,linewidth=2)
        plt.margins(x=0.1, y=0.1)
        plt.show()
        
        pass

    def tension_in_section1(Data,R1_tension):
        t1_data = []
        l1_data = []
        torque_data = []
        l1 = 100
        condition =  False 
        while condition == False:
            if ax.theta(Data,l1,R)*(180/np.pi) > 52:
                pass
            else:
                up_t1,*_ = ax.up_tension(Data,l1,R)
                t1,t2,t3,torque = ax.down_tension(Data,l1,R)
                t1_data.append(up_t1)
                l1_data.append(l1)
                torque_data.append(torque)

            if round(l1,1) == Data.max:
                condition = True
            else:
                l1 += 10

        custom_style = {
            'font.size': 16,  # Tamanho adequado para leitura de gráficos
            'axes.labelsize': 18,  # Tamanho dos rótulos dos eixos
            'axes.titlesize': 14,  # Tamanho do título do gráfico
            'axes.linewidth': 1.5,  # Espessura das bordas dos gráficos
            'xtick.labelsize': 18,  # Tamanho do texto dos ticks no eixo x
            'ytick.labelsize': 18,  # Tamanho do texto dos ticks no eixo y
            'lines.linewidth': 2,  # Espessura das linhas dos gráficos
            'lines.markersize': 6,  # Tamanho dos marcadores
            'legend.fontsize': 18,  # Tamanho da legenda
            'legend.frameon': False,  # Remove a moldura ao redor da legenda
            'legend.loc': 'best',  # Melhor posição automática para a legenda
            'figure.figsize': (4, 3),  # Tamanho padrão da figura (polegadas)
            'savefig.dpi': 1000,  # Alta resolução para exportação (publicação)
            'savefig.bbox': 'tight',  # Salva a imagem sem cortar parte do gráfico
            }


        plt.rcParams.update(custom_style)

        plt.plot(l1_data,t1_data, label=f'Tension with radius of ${R}m$', color="blue")
        plt.title("Tension behavior over the length of section one ")
        plt.xlabel("Length ($m$)")
        plt.ylabel("Tension ($N/m²$)")
        plt.legend()
        plt.grid(alpha=1,linewidth=2)
        plt.margins(x=0.1, y=0.1)
        plt.show()

        
        plt.plot(l1_data,torque_data, label=f'Torque with radius of ${R}m$', color="red", linestyle="--")
        plt.title("Torque behavior over the length of section one")
        plt.xlabel("Length ($m$)")
        plt.ylabel("Torque ($N*m$)")
        plt.legend()
        plt.grid(alpha=1,linewidth=2)
        plt.margins(x=0.1, y=0.1)
        plt.show()


        pass'''
