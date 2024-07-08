import csv
import matplotlib.pyplot as plt
displacement = []
force = []
with open('./data/force_disp_auxetic.csv', newline='') as csvfile:
    data_reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
    for row in data_reader:
        displacement.append(float(row[0][:10]))
        force.append(float(row[1]))

plt.figure(1)
plt.scatter(displacement, force)
plt.xlabel("Displacement [mm]")
plt.ylabel("Force [N]")
plt.grid()
plt.title("Auxetic Strip Force-Displacement Curve")
plt.savefig("example_strain_curve.png")