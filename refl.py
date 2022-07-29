import matplotlib.pyplot as plt
import numpy

df1 = open("watercorr_txt.txt").read().split("\n")
#df2 = open("reflectiont3conc19.txt").read().split("\n")
#df3 = open("reflectiont3conc31.txt").read().split("\n")
#df4 = open("reflectiont0conc01.txt").read().split("\n")

refl_list_1 = []
refl_list_2 = []
refl_list_3 = []
refl_list_4 = []
wl_list = []

for i in range(len(df1)):
    refl_list_1.append(df1[i])
#    refl_list_2.append(df2[i])
#    refl_list_3.append(df3[i])
#    refl_list_4.append(df4[i])

refl_list_1 = [float(e) for e in refl_list_1]
#refl_list_2 = [float(e) for e in refl_list_2]
#refl_list_3 = [float(e) for e in refl_list_3]
#refl_list_4 = [float(e) for e in refl_list_4]

for wl in range(380,781, 5):
    wl_list.append(wl)

plt.ylim(0.2, 1)
plt.plot(wl_list, refl_list_1, label='1% mel')
#plt.plot(wl_list, refl_list_2, '.', label='19% mel')
#plt.plot(wl_list, refl_list_3,'+', label='31% mel')
#plt.plot(wl_list, refl_list_4, '--', label='125 hrs')
plt.legend()
plt.ylabel('Reflectance R (a.u.)')
plt.xlabel('Wavelength (nm)')
plt.show()
