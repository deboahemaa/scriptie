import numpy as np
import matplotlib.pyplot as plt
import random
import math
import plotly.express as px
import pandas as pd

# parameters & lists ### DIT MOET IN FUNCTIE
#RI_list = [1, 1.5, 1.34, 1.40, 1.39, 1.40, 1.38, 1.44, 1] # Refractive index
#SC_list = [0, 100, 45, 30, 35, 25, 30, 5] # Scattering coefficient [mm-1]
#AC_list= [0, 0.00586, 0.00482, 0.03341, 0.24130, 0.03341, 0.08078, 0.04127] # Absorption coefficient [mm-1]
#anisotropy_list = [0, 0.86, 0.8, 0.9, 0.95, 0.8, 0.95, 0.75] # anisotropy [] 
#thickness_list = [0, 0.02, 0.11, 0.285, 0.375, 1.875, 1.975, 8.225] # mm
###### MOET IN FUNCTIE

test_ri = [1, 1, 1]
test_abs = [0, 10]
test_sc = [0, 90]
test_ani = [0, 0.75]
test_thick = [0, 0.02]

W_DR_list = [] # weight list for exiting photons
W_list_Abs = [] # weight list for photons at different depths
W_list_sc = [] 
W_list_T = []
W_list_R = []
layers = ["air", "Stratum corneum", "Living epidermis", "Papillary dermis", \
          "Upper blood net dermis", "reticular dermis", \
          "Deep blood net dermis", "subcutaneous fat"]

# Specular reflection when light propagates through one medium onto the tissue i.e. air to skin
def reflection_spec(n1, n2): 
    if n1 == n2:
        R_sp = 0
    else:
        R_sp = ((n1 - n2)**2) / ((n1+n2)**2) 
    return float(1 - R_sp) # Photon weight upon entering medium

def step_size(layer):
    i = layer
    ut = test_sc[i] + test_abs[i]
    ksi = random.uniform(0.00001 ,1)
    # Step size and probability distributions
    s_i = -np.log(ksi)/ut
    return s_i, ksi

# Function for propagation of single photon
def photon_moving(N):
    # Variables
    i = 1 #i-th layer of skin
    R_alpha = 0
    photon_count = 1
    dW = 0
    R_d = 0
    T_d = 0 
    W_start = reflection_spec(test_ri[0], test_ri[1])
    W = W_start


    # Variables and lists of coefficients changing per layer
    ut = test_sc[i] + test_abs[i] # [cm E-1] Scattering coefficients
    ua = test_abs[i] # [cm E-1] Absorption coefficients
    d_1 = test_thick[i]
    d_0 = test_thick[i-1]
    g = test_ani[i]
    d_b = d_1

    # Lists
    r_list = []
    z_list = [] # depth list
    coord = [0, 0, 0] # photon position in cartesian coordinates(x,y,z)
    mov = [0,0,0.999999] # directional cosines (ux, uy, uz)

    ## Photon propagation while weight is not too low 
    while W > 0.001 and photon_count <= N:

        # Step size determination process and determining concentrations of layer
        s_i = step_size(i, wl)  # save stepsize

        # Updated photon position     
        coord[0] = coord[0] + mov[0] * s_i/ut
        coord[1] = coord[1] + mov[1] * s_i/ut
        coord[2] = coord[2] + mov[2] * s_i/ut

        # lists for weight distribution graph 
        W_list_sc.append(W)
        r = np.sqrt(coord[0]**2 + coord[1]**2)
        r_list.append(r)
        z_list.append(coord[2])

        ### Prerequisites for transmittance and reflectance calcs
        # Statements to determine whether photon hits next layer during step
        if mov[2] < 0 and coord[2] >= test_thick[i-1]:
            d_b = (coord[2] - test_thick[i-1]) # distance to boundary
        elif mov[2] > 0 and coord[2] <= test_thick[i]:
            d_b = (test_thick[i] - coord[2]) # distance to boundary
        else:
            d_b = 0

        # Statements following the situation that particle might go to next layer
        if s_i >= d_b*ut:
            # Angle of incidince and calculation probability of internal reflection
            alpha_i = np.arccos(mov[2]) # angle of incidince
            
            if mov[2] >= 0:
                alpha_t = np.arcsin(np.sin(alpha_i) * test_ri[i] / test_ri[i+1])# tranmittance angle when going down
            elif mov[2] < 0:
                alpha_t = np.arcsin(np.sin(alpha_i) * test_ri[i] / test_ri[i-1])# tranmittance angle when going up
            
            # The angle of incidince greater than critical angle causes reflection so R = 1
            if mov[2] >= 0:
                if alpha_i > 0 and alpha_i < np.arcsin(test_ri[i+1]/test_ri[i]):
                    R_alpha = 0.5*((np.sin(alpha_i-alpha_t)/np.sin(alpha_t+alpha_i))**2 + (np.tan(alpha_i-alpha_t)/np.tan(alpha_t+alpha_i))**2)
                elif alpha_i == 0:
                    R_alpha = ((test_ri[i]-test_ri[i+1])*(test_ri[i]-test_ri[i+1]))/((test_ri[i]+test_ri[i+1])*(test_ri[i]+test_ri[i+1]))
                else:
                    R_alpha = 1

            elif mov[2] < 0: 

                if alpha_i > 0 and alpha_i < np.arcsin(test_ri[i-1]/test_ri[i]):
                    R_alpha = 0.5*((np.sin(alpha_i-alpha_t)/np.sin(alpha_t+alpha_i))**2 + (np.tan(alpha_i-alpha_t)/np.tan(alpha_t+alpha_i))**2)
                elif alpha_i == 0:
                    R_alpha = ((test_ri[i]-test_ri[i-1])*(test_ri[i]-test_ri[i-1]))/((test_ri[i]+test_ri[i-1])*(test_ri[i]+test_ri[i-1]))
                else:
                    R_alpha = 1
                
            ### Calculation of reflectance and transmittance
            # Reflectance at layer
            ksi_1 = ksi()
            if ksi_1 <= R_alpha and mov[2] >= 0:
                mov[2] = mov[2] * -1 
                coord[2] = test_thick[i]

            elif ksi_1 <= R_alpha and mov[2] < 0:
                mov[2] = mov[2] * -1 
                coord[2] = test_thick[i-1]

            # Succesful transmittance in next layer
            elif ksi_1 > R_alpha and mov[2] >= 0:
                coord[2] = test_thick[i] + (s_i - d_b*ut)
                
                i = i + 1
                
                if i <= len(test_thick)-1:    
                    # Variables and lists of coefficients changing per layer
                    ut = test_sc[i] + ua_layer(i, wl) # [mm E-1] Scattering coefficients
                    us = test_sc[i]
                    ua = ua_layer(i, wl) # [mm E-1] Absorption coefficients
                    g = test_ani[i]
                    
                    

                mov[0] = mov[0] * (np.sin(alpha_t)/np.sin(alpha_i))
                mov[1] = mov[1] * (np.sin(alpha_t)/np.sin(alpha_i))
                mov[2] = mov[2] * np.sin(alpha_t) / abs(mov[2])

            # Succesful transmittance to previous layer
            elif ksi_1 > R_alpha and mov[2] < 0:
                i = i - 1
                
                if i > 0:

                    # Variables and lists of coefficients changing per layer
                    ut = test_sc[i] + ua_layer(i, wl) # [mm E-1] Scattering coefficients
                    us = test_sc[i]
                    ua = ua_layer(i, wl) # [mm E-1] Absorption coefficients
                    g = test_ani[i]
                    coord[2] = test_thick[i] - (s_i - d_b*ut) 

                mov[0] = mov[0] * (np.sin(alpha_t)/np.sin(alpha_i))
                mov[1] = mov[1] * (np.sin(alpha_t)/np.sin(alpha_i))
                mov[2] = mov[2] * np.sin(alpha_t) / abs(mov[2])

            # Photon diffuse reflectance (Weight of photons exiting through top medium into detector)
            if i == 0 or coord[2]<=test_thick[0]:  
                R_d += W
                W_surf.append(W)
                x_list.append(coord[0])
                y_list.append(coord[1])
                W = W_start
                photon_count += 1
                step = 0
                coord = [0, 0, 0] # photon position in cartesian coordinates(x,y,z)
                mov = [0,0,0.999999] # directional cosines (ux, uy, uz)
                i = 1 # back to first layer 

                # Variables and lists of coefficients changing per layer
                ut = test_sc[i] + ua_layer(i, wl) # [mm E-1] Scattering coefficients
                us = test_sc[i]
                ua = ua_layer(i, wl) # [mm E-1] Absorption coefficients
                g = test_ani[i]
            
            # Photon diffuse transmittance(Weight of photons exiting through bottom medium)
            elif i > len(test_thick)-1 or coord[2] >= test_thick[-1]:
                T_d += W 
                W = W_start
                photon_count += 1
                step = 0
                coord = [0, 0, 0] # photon position in cartesian coordinates(x,y,z)
                mov = [0,0,0.999999] # directional cosines (ux, uy, uz)
                i = 1 # back to first layer

                # Variables and lists of coefficients changing per layer
                ut = test_sc[i] + ua_layer(i, wl) # [cm E-1] Scattering coefficients
                us = test_sc[i]
                ua = ua_layer(i, wl) # [cm E-1] Absorption coefficients
                g = test_ani[i]
            else:
                # New weight calculation if photon goes to next layer
                dW = W * (ua/(us+ua)) # ua absorption coefficient
                W_list_Abs.append(dW)
                A_d += dW
                W = W * (1 - ua/(us+ua))

        # if s< db, move photon and photon must experience absorption and scattering
        else:
            ### Photon scattering
            # angle calculation
            phi = 2*np.pi*ksi() # azimuthal angle phi [0,2pi]
            if g == 0:
                theta = np.arccos(2*ksi() - 1) # deflection angle theta [0,pi]
            else:
                theta = np.arccos((1/(2*g))*(1 + g**2 - ((1-g**2)/(1-g+2*g*ksi()))**2))

            # Calculating scatterdirection
            if abs(mov[2])>0.99999:
                mov[0] = np.sin(theta) * np.cos(phi)
                mov[1] = np.sin(theta) * np.sin(phi)
                mov[2] =  mov[2] * np.cos(theta) / abs(mov[2])

            else:
                mov[0] = np.sin(theta) * (mov[0]*mov[2]*np.cos(phi) - mov[1]*np.sin(phi)) / np.sqrt(1-mov[2]**2) + mov[0] * np.cos(theta)
                mov[1] = np.sin(theta) * (mov[1]*mov[2]*np.cos(phi) + mov[0]*np.sin(phi)) / np.sqrt(1-mov[2]**2) + mov[1] * np.cos(theta)
                mov[2] = -1* np.sin(theta) * np.cos(phi) * np.sqrt(1-mov[2]**2) + mov[2]*np.cos(theta)

            # New weight calculation if photon does not go to next layer
            dW = W * (ua/(us+ua)) # ua absorption coefficient
            W_list_Abs.append(dW)
            A_d += dW
            W = W * (1 - ua/(us+ua))

        # Photon diffuse transmittance(Weight of photons exiting through bottom medium)
        ### Russian roulette for photon 
        ksi_2 = ksi()
        if W<=0.001 and ksi_2 <=1/10:
            W = 10*W 

        # If photon dies, create new photon and return to start position with start values ## check photon afstand tot bron
        elif W<=0.001 and ksi_2 > 1/10:
            W = W_start
            photon_count += 1
            step = 0
            coord = [0, 0, 0] # photon position in cartesian coordinates(x,y,z)
            mov = [0,0,0.999999] # directional cosines (ux, uy, uz)
            i = 1

            # Variables and lists of coefficients changing per layer
            ut = test_sc[i] + ua_layer(i, wl) # [mm E-1] Scattering coefficients
            us = test_sc[i]
            ua = ua_layer(i, wl) # [mm E-1] Absorption coefficients
            g = test_ani[i]

    # Generate figure
    #print(R_d/N, T_d/N)
    #df = pd.DataFrame(np.column_stack((z_list, r_list)), columns=['Radius (cm)','Depth (cm)'])
    #fig = px.scatter(df, x="Radius (cm)", y="Depth (cm)", color=W_list_sc,
    #title="Numeric weight values mean continuous color")
    #fig.write_image("resultaat_test_mismatch.pdf")
    #plt.scatter(z_list, r_list, c= W_list_sc)
    #plt.show()
    return R_d/N, T_d/N


def Refl_spectrum(N):
    wl_start = 450
    wl = wl_start

    while wl <= 1050:
        coeff = Coeff_calc(wl)
        SC_list = coeff[0]
        AC_list = coeff[1]
        anisotropy_list = coeff[2]
        RI_list = coeff[3] 
        R_T_values = photon_moving(N)
        Ref = R_T_values[0]
        Tr = R_T_values[1]
        R_list.append(Ref)
        T_list.append(Tr)
        wl_list.append(wl)
        wl += 3
    plt.scatter(wl_list, R_list)
    plt.show()
    return

def average(k):
    Rd_avg = []
    Td_avg = []
    for i in range(k):
        x = photon_moving(50000)
        Rd_avg.append(x[0])
        Td_avg.append(x[1])
    print(np.mean(Rd_avg), np.std(Rd_avg), np.mean(Td_avg), np.std(Td_avg))
    return np.mean(Rd_avg), np.std(Rd_avg), np.mean(Td_avg), np.std(Td_avg)

average(10)