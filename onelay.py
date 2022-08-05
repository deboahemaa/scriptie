import numpy as np
import matplotlib.pyplot as plt
import random
import math
import plotly.express as px
import pandas as pd
import time
from colour_system import cs_srgb as cs
from PIL import Image
from PIL import ImageColor

### processing data to lists

# Data processing for water
wl_h2o_list = []
ua_h2o_list = []

df = open("data_water.csv").read().split("\n")
data = df[2:]

for i in range(len(data)):
    a, b = data[i].split("\t")
    wl_h2o_list.append(a)
    ua_h2o_list.append(b)

wl_h2o_list = [float(e) for e in wl_h2o_list]
ua_h2o_list = [float(e) for e in ua_h2o_list]

# Data processing for bilirubin
wl_bi_list = []
e_bi_list = []

df = open("data_bilirubin.csv").read().split("\n")
data = df[2:]

for i in range(len(data)):
    a, b = data[i].split("\t")
    wl_bi_list.append(a)
    e_bi_list.append(b)

wl_bi_list = [float(e) for e in wl_bi_list]
e_bi_list = [float(e) for e in e_bi_list]

# Data processing for deoxyhemoglobin and hemoglobin
e_hbo2_list = []
e_hb_list = []
wl_hb_list = []
arr = []
df = open("data_hemoglobin.csv").read().split("\n")
data = df[2:]

for i in range(len(data)-1):
    k = data[i].split("\t")
    arr.append(k)     
    wl_hb_list.append(arr[i][0])
    e_hbo2_list.append(arr[i][1])
    e_hb_list.append(arr[i][2])

wl_hb_list = [float(e) for e in wl_hb_list]
wl_hbo2_list = [float(e) for e in wl_hb_list]
e_hbo2_list = [float(e) for e in e_hbo2_list]
e_hb_list = [float(e) for e in e_hb_list]

# Lists of refractive index, scattering coefficient, anisotropy and thickness per layer
test_ri = [1.38, 1.44, 1] # 1, 1.5, 1.34, 1.4, 1.39, 1.4, 1.38, 1.44, 1
test_ani = [0.95, 0.75] # 0, 0.86, 0.8, 0.9, 0.95, 0.8, 0.95, 0.75
test_thick = [0, 0.6] # 0, 0.01, 0.07, 0.1, 0.1, 1.6, 0.12, 6

# Concetrations in different layers
# melanine, water, blood
def conc(i):
    if i == 2:
        C_list = [0, 0.05, 0]
    elif i == 3:
        C_list = [0.01, 0.2, 0]
    elif i == 4:
        C_list = [0, 0.5, 0.04]
    elif i == 5:
        C_list = [0, 0.6, 0.3]
    elif i ==6:
        C_list = [0, 0.7, 0.05]
    elif i ==7:
        C_list = [0, 0.7, 0.14]
    elif i ==1:
        C_list = [0, 0.65, 0.06]
    return C_list

# Relevant constants
S = 0.6
phi_hb = 0.25
H = 0.45
gamma = H*phi_hb
C_ext = 0.6
C_bi = 0

### Absorption coefficient for chromophores
#https://omlc.org/classroom/ece532/class3/muaspectra.html

#absorption coefficient of hemoglobin and waterfree tissue
def ua_other(wl):
    return 0.0244 + 8.53 / math.exp((wl-154)/66.2)

#absorption coefficient of melanin
def ua_me(wl):
    C_me, C_h2o, C_blood = conc(2)
    return ((2.37*10**4)*np.exp(-0.0056*wl) + (1.01*10**5)*np.exp(-0.0087*wl)) * np.log(10) * C_me

# Calculate absorption coefficient deoxyhemoglobin 
def e_hb(wl):
    if wl > 250:
        if wl in wl_hb_list:
            return e_hb_list[wl_hb_list.index(wl)]
        else:
            return np.mean(e_hb_list[wl_hb_list.index(wl+1)]+e_hb_list[wl_hb_list.index(wl-1)])

def ua_Hb(wl):
    if wl > 250:
        return 0.00054*e_hb(wl) 

# calculate absorption coefficient hemoglobin
def e_hbo2(wl):
    if wl > 250:
        if wl in wl_hbo2_list:
            return e_hbo2_list[wl_hbo2_list.index(wl)]
        else:
            return np.mean(e_hbo2_list[wl_hbo2_list.index(wl+1)]+e_hbo2_list[wl_hbo2_list.index(wl-1)])

def ua_hbo2(wl):
    if wl > 250:
        return 0.00054*e_hbo2(wl) 

# calculate absorption coefficient water
def ua_h2o(wl):
    if wl >= 225:
        if wl in wl_h2o_list:
            return ua_h2o_list[wl_h2o_list.index(wl)]
        else:
            wl = min(wl_h2o_list, key=lambda x:abs(x-wl))
            return np.mean(ua_h2o_list[wl_h2o_list.index(wl)] + ua_h2o_list[wl_h2o_list.index(wl) - 1])
       
# absorption coefficient for bilirubin
# molar concentration 0.05 to 0.4 g/L
def e_bi(wl):
    if wl in wl_bi_list:
        return e_bi_list[wl_bi_list.index(wl)]

def ua_bi(wl):
    if C_bi == 0:
        return 0
    if wl in wl_bi_list:
        return e_bi(wl) * C_bi * 0.00054

def test_absorp():
    C_me, C_h2o, C_blood = conc(2)

    ua_hbo2_list = []
    ua_hb_list = []
    ua_bi_list = []
    lambda_list = []
    ua_mel_list = []
    ua_other_list = []
    ua_h2o_list_test = []

    for wl in range(2, 1000, 2):
        lambda_list.append(wl)
        ua_hbo2_list.append(ua_hbo2(wl))
        ua_hb_list.append(ua_Hb(wl))
        ua_mel_list.append(ua_me(wl))
        ua_h2o_list_test.append(ua_h2o(wl))
        ua_bi_list.append(ua_bi(wl))
        ua_other_list.append(ua_other(wl))
    plt.plot(lambda_list, ua_bi_list, label='bi')
    plt.plot(lambda_list, ua_hbo2_list, label='hbo2')
    plt.plot(lambda_list, ua_hb_list, label='hb')
    plt.plot(lambda_list, ua_mel_list, label='mel')
    plt.plot(lambda_list, ua_other_list, label='other')
    plt.plot(lambda_list, ua_h2o_list_test, label='h2o')
    plt.yscale("log")
    plt.ylabel("absorption coefficient μ_a (mm E-1)")
    plt.xlabel("wavelength λ (nm)")
    plt.ylim(10**-5, 5*10**2)
    plt.legend()
    plt.show()
    return 

def test_refl():
    us_derm = []
    us_epi = []
    us_subcu = []
    lamb_list = []
    for wl in range(380, 780, 5):
        us_epi.append(us_epidermis(wl, 1))
        us_derm.append(us_dermis(wl, 3))
        us_subcu.append(us_subc(wl, 7))
        lamb_list.append(wl)
    plt.plot(lamb_list, us_epi, label='epidermis')
    plt.plot(lamb_list, us_derm, label='dermis')
    plt.plot(lamb_list, us_subcu, label='subcutaenous fat')
    plt.yscale("log")
    plt.ylabel("scattering coefficient μ_s (mm E-1)")
    plt.xlabel("wavelength λ (nm)")
    plt.legend()
    plt.show()
    return

#### Absorption coefficient for layers
#function to calculate absorption coefficient for stratum corneum
def ua_L1(wl, i):
    C_me, C_h2o, C_blood = conc(i)
    return ua_h2o(wl) * C_h2o + ua_other(wl)

# function to calculate absorption coefficient for living epidermis
def ua_L2(wl,i):
    C_me, C_h2o, C_blood = conc(i)
    return (C_me*ua_me(wl) + (1-C_me)*ua_other(wl))*(1-C_h2o) + C_h2o*ua_h2o(wl)

# function to calculate absorption coefficient for dermal layers
def ua_L3(wl, i):
    C_me, C_h2o, C_blood = conc(i)
    return (1-S)*gamma*(C_blood+C_ext)*ua_Hb(wl) + ua_hbo2(wl)*gamma*(C_blood+C_ext)*S + \
ua_bi(wl)*C_bi*(1-gamma*(C_blood+C_ext))+ ua_h2o(wl) *C_h2o*(1-gamma*(C_blood+C_ext))*(1-C_bi) + \
ua_other(wl)*(1-gamma*(C_blood+C_ext))*(1-C_h2o)*(1-C_bi)

# Define total absorption coefficient based on layer
def ua_layer(i, wl):
    if i == 3:
        ua = ua_L1(wl, i)
    elif i == 2:
        ua = ua_L2(wl, i)
    elif i==1:
        ua = ua_L3(wl, i)
    return ua

### Scattering coefficient for layers
# function to calculate scattering coefficient for epidermal layers
def us_epidermis(wl, i):
    return 68.7 * (wl/500)**(-1.161)

# function to calculate scattering coefficient for dermal layers
def us_dermis(wl, i):
    return 45.3*(wl/500)**(-1.292)

# function to calculate scattering coefficient for subcutaneous fat
def us_subc(wl, i):
    return 35.2*(wl/500)**(-0.988)

# Define total absorption coefficient based on layer
def us_layer(i, wl):
    if i == 3 or i == 2:
        us = us_epidermis(wl, i)
    elif i == 7:
        us = us_dermis(wl, i)
    elif i==1:
        us = us_subc(wl, i)
    return us

### The monte carlo simulations
W_list_Abs = [] # weight list for photons at different depths
W_list_sc = []
W_surf = []
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

# Random number uniform between 0 and 1
def ksi():
    return random.uniform(0,1)

# Step size function
def step_size(layer, wl):
    i = layer
    ut = us_layer(i, wl) + ua_layer(i, wl)
    s_i = -np.log(ksi())/ut
    return s_i

# Function for propagation of single photon
def photon_moving(N, wl):
    # Variables
    i = 1 #i-th layer of skin
    R_alpha = 1
    photon_count = 1
    dW = 0
    R_d = 0
    T_d = 0
    A_d = 0 
    W_start = float(reflection_spec(test_ri[0], test_ri[1]))
    W = 1

    # Variables and lists of coefficients changing per layer
    us = us_layer(i, wl)
    g = test_ani[i]
    d_b = test_thick[i]
    C_me, C_h2o, C_blood = conc(i)
    ua = ua_layer(i, wl) # [cm E-1] Absorption coefficients
    ut = us+ua # [cm E-1] Scattering coefficients

    # Lists
    x_list = []
    y_list = []
    r_list = []
    z_list = [] # depth list
    coord = [0, 0, 0] # photon position in cartesian coordinates(x,y,z)
    mov = [0,0,0.999999] # directional cosines (ux, uy, uz)

    ## Photon propagation while weight is not too low 
    while W > 0.001 and photon_count <= N:

        # Step size determination process and determining concentrations of layer
        s_i = step_size(i, wl)  # save stepsize
        C_me, C_h2o, C_blood = conc(i)

        # Updated photon position    
        coord[0] = coord[0] + mov[0] * s_i
        coord[1] = coord[1] + mov[1] * s_i
        coord[2] = coord[2] + mov[2] * s_i/ut

#        if abs(coord[0]) <= 0.1 and abs(coord[1]) <=0.1:
#            C_ext = 0.1
#            S = 0.4
#        else:
#            S = 0.4
#            C_ext = 0

#        if abs(coord[0]) <= 0.20 and abs(coord[1]) <=0.2:
#            C_bi = 0.35
#            S = 0.4
#        else:
#            C_bi = 0
#            S = 0.4
        
        # lists for weight distribution graph  
        W_list_sc.append(W)
        r_list.append(coord[1])
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
                if alpha_i == 0:
                    R_alpha = ((test_ri[i]-test_ri[i+1])*(test_ri[i]-test_ri[i+1]))/((test_ri[i]+test_ri[i+1])*(test_ri[i]+test_ri[i+1]))
                elif alpha_i-alpha_t == 0:
                    R_alpha = 0
                elif alpha_i > 0 and alpha_i < np.arcsin(test_ri[i+1]/test_ri[i]):
                    R_alpha = 0.5*((np.sin(alpha_i-alpha_t)/np.sin(alpha_t+alpha_i))**2 + (np.tan(alpha_i-alpha_t)/np.tan(alpha_t+alpha_i))**2)
                else:
                    R_alpha = 1

            elif mov[2] < 0:
                if alpha_i == 0:
                    R_alpha = ((test_ri[i]-test_ri[i-1])*(test_ri[i]-test_ri[i-1]))/((test_ri[i]+test_ri[i-1])*(test_ri[i]+test_ri[i-1])) 
                elif alpha_i-alpha_t == 0:
                    R_alpha = 0
                elif alpha_i > 0 and alpha_i < np.arcsin(test_ri[i-1]/test_ri[i]):
                    R_alpha = 0.5*((np.sin(alpha_i-alpha_t)/np.sin(alpha_t+alpha_i))**2 + (np.tan(alpha_i-alpha_t)/np.tan(alpha_t+alpha_i))**2)
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
                    C_me, C_h2o, C_blood = conc(i)
                    us = us_layer(i, wl)
                    ua = ua_layer(i, wl) # [mm E-1] Absorption coefficients
                    ut = us+ua # [mm E-1] Scattering coefficients
                    g = test_ani[i] 

                mov[0] = mov[0] * (np.sin(alpha_t)/np.sin(alpha_i))
                mov[1] = mov[1] * (np.sin(alpha_t)/np.sin(alpha_i))
                mov[2] = mov[2] * np.sin(alpha_t) / abs(mov[2])

            # Succesful transmittance to previous layer
            elif ksi_1 > R_alpha and mov[2] < 0:
                i = i - 1
                
                if i > 0:

                    # Variables and lists of coefficients changing per layer
                    C_me, C_h2o, C_blood = conc(i)
                    us = us_layer(i, wl)
                    ua = ua_layer(i, wl) # [mm E-1] Absorption coefficients
                    ut = us+ua # [mm E-1] Scattering coefficients
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
                W = 1
                photon_count += 1
                step = 0
                coord = [0, 0, 0] # photon position in cartesian coordinates(x,y,z)
                mov = [0,0,0.999999] # directional cosines (ux, uy, uz)
                i = 1 # back to first layer 

                # Variables and lists of coefficients changing per layer
                us = us_layer(i, wl)
                ua = ua_layer(i, wl) # [mm E-1] Absorption coefficients
                ut = us+ua # [mm E-1] Scattering coefficients
                g = test_ani[i]
            
            # Photon diffuse transmittance(Weight of photons exiting through bottom medium)
            elif i > len(test_thick)-1 or coord[2] >= test_thick[-1]:
                T_d += W 
                W = 1
                photon_count += 1
                step = 0
                coord = [0, 0, 0] # photon position in cartesian coordinates(x,y,z)
                mov = [0,0,0.999999] # directional cosines (ux, uy, uz)
                i = 1 # back to first layer

                # Variables and lists of coefficients changing per layer
                us = us_layer(i, wl)
                ua = ua_layer(i, wl) # [cm E-1] Absorption coefficients
                ut = us+ua # [mm E-1] Scattering coefficients
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
            W = 1
            photon_count += 1
            step = 0
            coord = [0, 0, 0] # photon position in cartesian coordinates(x,y,z)
            mov = [0,0,0.999999] # directional cosines (ux, uy, uz)
            i = 1

            # Variables and lists of coefficients changing per layer
            us = us_layer(i, wl)
            ua = ua_layer(i, wl) # [mm E-1] Absorption coefficients
            ut = us+ua # [mm E-1] Scattering coefficients
            g = test_ani[i]

    # Generate values
    print(R_d/N, T_d/N)
    #plt.scatter(r_list, z_list, c= W_list_sc)
    #plt.colorbar()
    #plt.show()
    return R_d/N, T_d/N, A_d/N, x_list, y_list, W_surf

def refl_spec(wl_start, wl_end, N):
    W_list_T = []
    W_list_R = []
    wl_list_it = []
    spec_img = []
    x_pos = 0
    y_pos = 0

    for wavelength in np.arange(wl_start, wl_end+1, 5):
        wl = wavelength
        R_d, T_d, A_d, x, y, W = photon_moving(N, wl)
        W_list_R.append(R_d)
        W_list_T.append(T_d)
        wl_list_it.append(wl)

        grid = np.zeros((25,25), float)
        counter = np.zeros((25,25))
        
        for i in range(len(x)):
            if abs(x[i])<=0.5 and abs(y[i])<=0.5: 
                x_pos = int((x[i] + 0.5)*25) 
                y_pos = int((y[i] + 0.5)*25)
                grid[x_pos, -y_pos] += W[i]
                counter[x_pos, -y_pos] += 1

        for x in range(25):
            for y in range(25):
                if counter[x,y] != 0:
                    grid[x,y] = grid[x,y] / counter[x,y]
        spec_img.append(grid)
    np.savetxt('reflection_L7', W_list_R)
    np.savetxt('transm_L7', W_list_T)
    plt.plot(wl_list_it, W_list_R)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Reflectance (A.U.)')
    plt.title('Reflectance spectrum of Subcutaneous Fat')
    plt.show()
    return W_list_T, W_list_R, wl_list_it, spec_img

def image():
    trans, refl, golf, spectral = refl_spec(380, 780, 1000) # wl_begin, wl_end, N_photons
    img = Image.new('RGB', (25, 25)) # 25 x 25 picture

    #np.savetxt("refl_spect_pixel_t0conc1", spectral) # save spectral images
    np.savetxt("melcorr1", refl) # save reflectance per wl in range (wl_begin, wl_end, 5)

    # generate picture from spectral images
    for x in range(25):
        for y in range(25):
            spec = []
            for lam in np.arange(380., 781., 5):
                spec_img = spectral[golf.index(lam)]
                spec.append(spec_img[x,y])

            color = cs.spec_to_rgb(np.asarray(spec), out_fmt='html')
            img.putpixel((x,y), ImageColor.getrgb(color))

    img.save('melcorr1.png')
    img.show() 
    return img

# Measure time it takes to run full cycle by running below: 
start = time.time()
#photon_moving(5000, 400)
#test_absorp()
#test_refl()
#image()
refl_spec(380,780,1000)
end = time.time()
print(end - start)