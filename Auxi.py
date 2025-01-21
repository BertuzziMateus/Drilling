import numpy as np
import sympy as sp


def signal_change(s)->tuple: #return arrays with the signals of the value, -1 for negative and 1 for positive
    a = None
    codition = None
    cs = None
    for i in range(len(s)):
        if i == 0 :
            if s[i] > 0 :
                codition = True
            else: 
                codition = False
        else:
            if s[i] > 0 :
                a = True
                if codition != a:
                    value = 'Yes'
                    cs = i
                    break
                else:
                    value = 'No'
            else: 
                a = False
                if codition != a:
                    value = 'Yes'
                    cs = i
                    break
                else:
                    value = 'No'
    return (value, cs)

def theta( Data, l1, R ) ->  float: # Estimated final angle in the curve
    l1 = l1
    x = sp.Symbol('x')
    x3 , y3 = Data.P3
    l3 = np.sqrt( ((x3- R)**2) + (y3-l1)**2 - (R**2) )
    a1 = (2*np.arctan((R - np.sqrt(R**2 - l1**2 + 2*l1*y3 + l3**2 - y3**2))/(-l1 + l3 + y3)))
    a2 = (2*np.arctan((R + np.sqrt(R**2 - l1**2 + 2*l1*y3 + l3**2 - y3**2))/(-l1 + l3 + y3)))
    if R == x3:
         a = a1
    else:
        f1 = ( -R*sp.cos(x) + l3*sp.sin(x) )
        f2 = ( l3*sp.cos(x) + R*sp.sin(x) )
        if (np.sign((f1.evalf(subs={x: a1})))) == (np.sign((x3-R))) and (np.sign((f2.evalf(subs={x: a1})))) == (np.sign(((y3-l1)))):
            a = a1
        else:
            a = a2
    return float(a)

def lenght( Data, l1,R )  -> list: # Returns l1,l2 and l3 lengths from the each part

    l3 = np.sqrt( ((Data.P3[0]-R)**2) + (Data.P3[1]-l1)**2 - (R**2) )

    angle  = theta(Data, l1,R)

    l2 = (2*np.pi*R) * (((angle)*(180/np.pi))/360)

    return l1,l2,l3

def curve_points( Data,l1,R ) -> list:  # Estimated points in the curve part
    angle = np.linspace(np.pi, np.pi + theta(Data, l1, R),1000)
    l1, *_ = lenght(Data,l1,R)
    y = l1 -  R * np.sin(angle)
    x = R * (np.cos(angle)) + R
    return x,y

def points_coordinates( Data,l1,R) -> list : # Estimated points in well
    x,y = [Data.P0[0], 0],[Data.P0[1],l1]

    points = curve_points(Data,l1,R)

    for i in range(len(points[0])):
            x.append(points[0][i])
            y.append(points[1][i])
    x.append(Data.P3[0]),y.append(Data.P3[1])
    
    return x,y 

def buckling( Data, l1, R ) -> float:

    PSB = (5000*4.4482216172334)*(Data.d_ext_command) # LBF - Newton

    alpha =  1 - ((Data.ro_fluid) / (Data.ro_command))
    ws = Data.lambd_command

    lc = round((PSB*1.2) / (ws*alpha*np.cos(theta(Data, l1, R))))
    
    codition =  True

    while codition == True:

        if (lc % 9) >= 5 and (lc % 9) != 0:
            # print(lc,'dif')
            # print(lc%9)
            # # if lc
            lc += 1
        elif (lc % 9) < 5 and (lc % 9) != 0:
            # print(lc)
            # print(lc%9)
            lc -= 1
        else:
            # print(lc)
            codition = False
    
    #print(lc,l1,R)

    return lc
    
def Nl( Data, l1, R ) -> float:
    lc = buckling(Data,l1,R)
    angle = theta(Data, l1, R)
    Y_F =  Data.P3[1]
    ff = ( Data.ro_fluid*Data.g*Y_F*(np.pi*(Data.d_ext_command**2 - Data.d_int_command**2  ) / 4))

    neutral_line = (Data.z + ff ) / ( Data.lambd_command* Data.g* ( np.cos(angle) - Data.µ* (1 - (Data.ro_fluid/Data.ro_command) )* np.sin(angle)  ) )

    if neutral_line < lc: # Neutral line is in command
   
        neutral_line = neutral_line
    
    elif neutral_line > lc: # Neutral line is in heavy pipe 
        Y_F =  Data.P3[1]
        Y_ic = Y_F - ( lc*np.cos(angle) )
        fic = ( Data.ro_fluid*Data.g*Y_ic*np.pi*( Data.d_ext_command**2 - Data.d_ext_heavy**2 ) ) / 4

        A = ( ff + Data.z - fic )  / ( Data.lambd_heavy* Data.g* ( np.cos(angle) - Data.µ* (1 - (Data.ro_fluid/Data.ro_heavypipe) )* np.sin(angle)  ) ) # frist term
        B = lc* ( (Data.lambd_command - Data.lambd_heavy)*np.cos(angle) - (Data.µ*( (Data.lambd_command*(1 - (Data.ro_fluid/Data.ro_command))) - (Data.lambd_heavy*(1 - (Data.ro_fluid/Data.ro_heavypipe))) ))*np.sin(angle) ) # second term
        C = Data.lambd_heavy*( np.cos(angle) - Data.µ* (1 - (Data.ro_fluid/Data.ro_heavypipe) )* np.sin(angle)  ) #third term
    

        neutral_line = A - (B/C)
    
        if neutral_line > (lc + Data.lp): # Neutral line is in drill
            Y_F =  Data.P3[1]
            Y_ic = Y_F - ( lc*np.cos(angle) )
            Y_ip = Y_F - ( ( lc + Data.lp )*np.cos(angle) )
            fic = ( Data.ro_fluid*Data.g*Y_ic*np.pi*( Data.d_ext_command**2 - Data.d_ext_heavy**2 ) ) / 4
            fip = ( Data.ro_fluid*Data.g*Y_ip*np.pi*( Data.d_ext_heavy**2 - Data.d_ext_drill**2 ) ) / 4

            A = ( ff + Data.z - fic  - fip )  / ( Data.lambd_drill* Data.g* ( np.cos(angle) - Data.µ* (1 - (Data.ro_fluid/Data.ro_drillpipe) )* np.sin(angle)  ) ) # frist term
            B = lc* ( (Data.lambd_command - Data.lambd_drill)*np.cos(angle) - (Data.µ*( (Data.lambd_command*(1 - (Data.ro_fluid/Data.ro_command))) - (Data.lambd_drill*(1 - (Data.ro_fluid/Data.ro_drillpipe))) ))*np.sin(angle) ) #second term
            C = Data.lambd_drill*( np.cos(angle) - Data.µ* (1 - (Data.ro_fluid/Data.ro_drillpipe) )* np.sin(angle)  ) #third term
            D = Data.lp* ( (Data.lambd_heavy - Data.lambd_drill)*np.cos(angle) - (Data.µ*( (Data.lambd_heavy*(1 - (Data.ro_fluid/Data.ro_heavypipe))) - (Data.lambd_drill*(1 - (Data.ro_fluid/Data.ro_drillpipe))) ))*np.sin(angle) ) #fourth term
         
            neutral_line = A - (B/C) - (D/C)

    return neutral_line

def up_tension( Data, l1, R) -> list : # Estimated tensions in upping

    lc = buckling(Data,l1,R)
    l1,l2,l3   =  lenght(Data, l1, R)
    angle  =  theta(Data, l1, R)
    weight_3 = ( ( (Data.lambd_command*lc) + (Data.lambd_heavy*Data.lp) + (Data.lambd_drill*(l3 - lc - Data.lp)) )*Data.g )
    buoyancy_3 = (( (Data.area_drill*(l3-Data.lp-lc)*Data.ro_fluid) + (Data.area_command*lc*Data.ro_fluid) + (Data.area_heavy*Data.lp*Data.ro_fluid) ) * Data.g)     
    

    Y_F =  Data.P3[1]
    Y_ic = Y_F - ( lc*np.cos(angle) )
    Y_ip = Y_F - ( ( lc + Data.lp )*np.cos(angle) )
    fic = ( Data.ro_fluid*Data.g*Y_ic*np.pi*( Data.d_ext_command**2 - Data.d_ext_heavy**2 ) ) / 4
    fip = ( Data.ro_fluid*Data.g*Y_ip*np.pi*( Data.d_ext_heavy**2 - Data.d_ext_drill**2 ) ) / 4
    ff = ( Data.ro_fluid*Data.g*Y_F*(np.pi*(Data.d_ext_command**2 - Data.d_int_command**2  ) / 4))


    tension_3 = weight_3*np.cos(angle) +  Data.µ*(weight_3 - buoyancy_3)*np.sin(angle) + fic + fip - ff

    angle_variation = np.arange(0.001,angle,0.001)[::-1]


    DNS = []

    for i in range(len( angle_variation)):

        DN = tension_3*angle_variation[i] - ( 1 - (Data.ro_fluid/Data.ro_drillpipe) )*Data.g*Data.lambd_drill*R*np.sin(angle_variation[i])*angle_variation[i]
        DNS.append(DN)
    
    signals = np.sign(DNS)


    coditions = signal_change(signals)

    condition_signal_change = coditions[0]
    condition_where_change = coditions[1]
 
 
    if condition_signal_change == "No":
        
        if DNS[0] > 0: # DN POSTIVO SEM TROCA DE SINAL
 
            A = Data.lambd_drill*Data.g*R
            C1 = ( A / ( (1 + Data.µ**2) ) ) * ( ( (Data.µ**2) * ( 1 - (Data.ro_fluid/Data.ro_drillpipe) ) )  - 1 )
            C2 =  - ( (A*Data.µ) / ( (1 + Data.µ**2) ) ) * ( 2 - (Data.ro_fluid/Data.ro_drillpipe) )
            
            K = ( tension_3 - C1*np.sin(angle) - C2*np.cos(angle) )  / ( np.exp(-Data.µ*angle) )


            tension_2 = C1*np.sin(0) + C2*np.cos(0) + K*( np.exp(-Data.µ*0) )


        else: # DN NEGATIVO SEM TROCA DE SINAL

            A = Data.lambd_drill*Data.g*R
            B1 = ( A / ( (1 + Data.µ**2) ) ) * ( ( (Data.µ**2) * ( 1 - (Data.ro_fluid/Data.ro_drillpipe) ) )  - 1 )
            B2 =  ( (A*Data.µ) / ( (1 + Data.µ**2) ) ) * ( 2 - (Data.ro_fluid/Data.ro_drillpipe) )
            K = ( tension_3 - B1*np.sin(angle) - B2*np.cos(angle) )  / ( np.exp(Data.µ*angle) )


            tension_2 = B1*np.sin(0) + B2*np.cos(0) + K*( np.exp(Data.µ*0) )


    if condition_signal_change == "Yes":


        A = Data.lambd_drill*Data.g*R
        B1 = ( A / ( (1 + Data.µ**2) ) ) * ( ( (Data.µ**2) * ( 1 - (Data.ro_fluid/Data.ro_drillpipe) ) )  - 1 )
        B2 =  ( (A*Data.µ) / ( (1 + Data.µ**2) ) ) * ( 2 - (Data.ro_fluid/Data.ro_drillpipe) )
        K = ( tension_3 - B1*np.sin(angle) - B2*np.cos(angle) )  / ( np.exp(Data.µ*angle) )


        tension_change = B1*np.sin(angle_variation[condition_where_change]) + B2*np.cos(angle_variation[condition_where_change]) + K*( np.exp(Data.µ*angle_variation[condition_where_change]) )

        C1 = ( A / ( (1 + Data.µ**2) ) ) * ( ( (Data.µ**2) * ( 1 - (Data.ro_fluid/Data.ro_drillpipe) ) )  - 1 )
        C2 =  - ( (A*Data.µ) / ( (1 + Data.µ**2) ) ) * ( 2 - (Data.ro_fluid/Data.ro_drillpipe) )
        
        K = ( tension_change - C1*np.sin(angle) - C2*np.cos(angle) )  / ( np.exp(-Data.µ*angle) )


        tension_2 = C1*np.sin(0) + C2*np.cos(0) + K*( np.exp(-Data.µ*0) )

   
    tension_1 = tension_2 + ( Data.lambd_drill*l1 )*Data.g

    
    return tension_1, tension_2, tension_3

def down_tension( Data, l1, R ) -> list : # Estimated tensions in descent
    lc = buckling(Data,l1,R)
    l1,l2 ,l3   =  lenght(Data, l1,R)
    angle  =  theta(Data, l1,R)
    ld = l3 - Data.lp - lc
    weight_3 = ( ( (Data.lambd_command*lc) + (Data.lambd_heavy*Data.lp) + (Data.lambd_drill*(l3 - lc - Data.lp)) )*Data.g )

    buoyancy_3 = (
        ( (Data.area_drill*(l3-Data.lp-lc)*Data.ro_fluid) + 
        (Data.area_command*lc*Data.ro_fluid) + 
        (Data.area_heavy*Data.lp*Data.ro_fluid) ) *
        Data.g
        )      
    
    Y_F =  Data.P3[1]
    Y_ic = Y_F - ( lc*np.cos(angle) )
    Y_ip = Y_F - ( ( lc + Data.lp )*np.cos(angle) )

    ff = ( Data.ro_fluid*Data.g*Y_F*(np.pi*(Data.d_ext_command**2 - Data.d_int_command**2  ) / 4))
    fic = ( Data.ro_fluid*Data.g*Y_ic*np.pi*( Data.d_ext_command**2 - Data.d_ext_heavy**2 ) ) / 4
    fip = ( Data.ro_fluid*Data.g*Y_ip*np.pi*( Data.d_ext_heavy**2 - Data.d_ext_drill**2 ) ) / 4

    tension_3 = weight_3*(  np.cos(angle) - Data.µ*np.sin(angle)  ) + ( Data.µ*buoyancy_3*np.sin(angle) ) + fic + fip - Data.z - ff

    angle_variation = np.arange(0.01,angle,0.01)[::-1]


    DNS = []

    for i in range(len( angle_variation)):

        DN = tension_3*angle_variation[i] - ( 1 - (Data.ro_fluid/Data.ro_drillpipe) )*Data.g*Data.lambd_drill*R*np.sin(angle_variation[i])*angle_variation[i]
        if DN < 0:
            a = 1
        DNS.append(DN)
    
    signals = np.sign(DNS)


    coditions = signal_change(signals)

    condition_signal_change = coditions[0]
    condition_where_change = coditions[1]
 

    if condition_signal_change == "No":
        
        if DNS[0] > 0: # DN POSTIVO SEM TROCA DE SINAL
           
            A = Data.lambd_drill*Data.g*R

            B1 = ( A / ( (1 + Data.µ**2) ) ) * ( ( (Data.µ**2) * ( 1 - (Data.ro_fluid/Data.ro_drillpipe) ) )  - 1 )  # PG - 4* 
            B2 = ( (A*Data.µ) / ( (1 + Data.µ**2) ) ) * ( 2 - (Data.ro_fluid/Data.ro_drillpipe) ) # PG - 4*
                        
            K = ( tension_3 - B1*np.sin(angle) - B2*np.cos(angle) )  / ( np.exp(Data.µ*angle) )

            tension_2 = B1*np.sin(0) + B2*np.cos(0) + K*( np.exp(Data.µ*0) )

            B3 = B1 - ( (1 - (Data.ro_fluid/Data.ro_drillpipe) )* Data.lambd_drill* Data.g*R )

            medial_N_force = B3*(1-np.cos(angle)) + B2*np.sin(angle) + (K/Data.µ)*( (np.exp(Data.µ*angle)) - 1 )


            frictional_torque_2 = ( Data.µ * medial_N_force* Data.d_ext_drill ) / 2

            fat3_command = Data.µ* ( 1 - (Data.ro_fluid/Data.ro_command) )* Data.lambd_command* lc*Data.g* np.sin(angle)
            fat3_heavy = Data.µ* ( 1 - (Data.ro_fluid/Data.ro_heavypipe) )* Data.lambd_heavy* Data.lp*Data.g* np.sin(angle)
            fat3_drill = Data.µ* ( 1 - (Data.ro_fluid/Data.ro_drillpipe) )* Data.lambd_drill* ld*Data.g* np.sin(angle)


            torque3_command = (Data.d_ext_command/2) * fat3_command
            torque3_heavy = (Data.d_ext_heavy/2) * fat3_heavy
            torque3_drill = (Data.d_ext_drill/2) * fat3_drill

            frictional_torque_3 = torque3_command + torque3_heavy  + torque3_drill

            torque = frictional_torque_2 + frictional_torque_3

        else: # DN NEGATIVO SEM TROCA DE SINAL

            A = Data.lambd_drill*Data.g*R
            C1 = ( A / ( (1 + Data.µ**2) ) ) * ( ( (Data.µ**2) * ( 1 - (Data.ro_fluid/Data.ro_drillpipe) ) )  - 1 )
            C2 =  - ( (A*Data.µ) / ( (1 + Data.µ**2) ) ) * ( 2 - (Data.ro_fluid/Data.ro_drillpipe) )
                    
            K = ( tension_3 - C1*np.sin(angle) - C2*np.cos(angle) )  / ( np.exp(-Data.µ*angle) )

            tension_2 = C1*np.sin(0) + C2*np.cos(0) + K*( np.exp(-Data.µ*0) )

            C3 = C1 - ( (1 - (Data.ro_fluid/Data.ro_drillpipe) )* Data.lambd_drill* Data.g*R )

            medial_N_force = C3*(1-np.cos(angle)) + C2*np.sin(angle) - (K/Data.µ)*( (np.exp(-Data.µ*angle)) - 1 )


            frictional_torque_2 = -(( Data.µ * medial_N_force* Data.d_ext_drill ) / 2)

                        
            fat3_command = Data.µ* ( 1 - (Data.ro_fluid/Data.ro_command) )* Data.lambd_command* lc*Data.g* np.sin(angle)
            fat3_heavy = Data.µ* ( 1 - (Data.ro_fluid/Data.ro_heavypipe) )* Data.lambd_heavy* Data.lp*Data.g* np.sin(angle)
            fat3_drill = Data.µ* ( 1 - (Data.ro_fluid/Data.ro_drillpipe) )* Data.lambd_drill* ld*Data.g* np.sin(angle)


            torque3_command = (Data.d_ext_command/2) * fat3_command
            torque3_heavy = (Data.d_ext_heavy/2) * fat3_heavy
            torque3_drill = (Data.d_ext_drill/2) * fat3_drill

            frictional_torque_3 = torque3_command + torque3_heavy  + torque3_drill

            torque = frictional_torque_2 + frictional_torque_3


    elif condition_signal_change == "Yes": # negative - positive

            
        A = Data.lambd_drill*Data.g*R
        C1 = ( A / ( (1 + Data.µ**2) ) ) * ( ( (Data.µ**2) * ( 1 - (Data.ro_fluid/Data.ro_drillpipe) ) )  - 1 )
        C2 =  - ( (A*Data.µ) / ( (1 + Data.µ**2) ) ) * ( 2 - (Data.ro_fluid/Data.ro_drillpipe) )
                
        K = ( tension_3 - C1*np.sin(angle) - C2*np.cos(angle) )  / ( np.exp(-Data.µ*angle) )

        C3 = C1 - ( (1 - (Data.ro_fluid/Data.ro_drillpipe) )* Data.lambd_drill* Data.g*R )

        medial_N_force  = -C3*( np.cos(angle_variation[condition_where_change])-np.cos(angle) ) + C2*(np.sin(angle_variation[condition_where_change]) - np.sin(angle)) - (K/Data.µ)*(np.exp(-Data.µ*angle_variation[condition_where_change])-np.exp(-Data.µ*angle))

        
        frictional_torque_2_change = - (( Data.µ * medial_N_force* Data.d_ext_drill ) / 2)


        # print(angle*57.29,'angulo')
        # print(angle_variation[condition_where_change]*57.29,'angulo de variação')


        tension_change = C1*np.sin(angle_variation[condition_where_change]) + C2*np.cos(angle_variation[condition_where_change]) + K*( np.exp(-Data.µ*(angle_variation[condition_where_change])) ) #CALCULA A TENSÃO NO PONTO QUE TROCOU O SINAL

    
        B1 = ( A / ( (1 + Data.µ**2) ) ) * ( ( (Data.µ**2) * ( 1 - (Data.ro_fluid/Data.ro_drillpipe) ) )  - 1 )
        B2 =  ( (A*Data.µ) / ( (1 + Data.µ**2) ) ) * ( 2 - (Data.ro_fluid/Data.ro_drillpipe) )
        K = ( tension_change - B1*np.sin(angle_variation[condition_where_change]) - B2*np.cos(angle_variation[condition_where_change]) )  / ( np.exp(Data.µ*angle_variation[condition_where_change]) )

        B3 = B1 - ( (1 - (Data.ro_fluid/Data.ro_drillpipe) )* Data.lambd_drill* Data.g*R )

        medial_N_force = B3*(1-np.cos(angle_variation[condition_where_change])) + B2*np.sin(angle_variation[condition_where_change]) + (K/Data.µ)*( (np.exp(Data.µ*(angle_variation[condition_where_change]))) - 1 )


        frictional_torque_2 = (( Data.µ * medial_N_force* Data.d_ext_drill ) / 2) + frictional_torque_2_change

    
                    
        fat3_command = Data.µ* ( 1 - (Data.ro_fluid/Data.ro_command) )* Data.lambd_command* lc*Data.g* np.sin(angle)
        fat3_heavy = Data.µ* ( 1 - (Data.ro_fluid/Data.ro_heavypipe) )* Data.lambd_heavy* Data.lp*Data.g* np.sin(angle)
        fat3_drill = Data.µ* ( 1 - (Data.ro_fluid/Data.ro_drillpipe) )* Data.lambd_drill* ld*Data.g* np.sin(angle)


        torque3_command = (Data.d_ext_command/2) * fat3_command
        torque3_heavy = (Data.d_ext_heavy/2) * fat3_heavy
        torque3_drill = (Data.d_ext_drill/2) * fat3_drill

        frictional_torque_3 = torque3_command + torque3_heavy  + torque3_drill

        torque = frictional_torque_2 + frictional_torque_3


        tension_2 = B1*np.sin(0) + B2*np.cos(0) + K*( np.exp(Data.µ*0) )






   
    tension_1 = tension_2 + ( Data.lambd_drill*l1 )*Data.g

    return tension_1, tension_2, tension_3,torque#


