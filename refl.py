import matplotlib.pyplot as plt
import numpy

df = open("raw_data_refl_def.csv").read().split("\n")

refl_list = []
wl_list = []

print(df[1].split(" "))
for i in range(len(df)):
    a, b = df[i].split(" ")
    refl_list.append(a)


refl_list = [float(e) for e in refl_list]
print(refl_list)

for wl in range(455,700, 5):
    wl_list.append(wl)

plt.ylim(0.2, 1)
plt.plot(wl_list, refl_list)
plt.show()