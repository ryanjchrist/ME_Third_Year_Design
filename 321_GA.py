import math
import numpy as np
from sympy import symbols, Eq, solve
import matplotlib.pyplot as plt
import pandas as pd
import itertools

def calculate_fluid_mechanical_power(Q_ft, rho, g, h, n_p):
    P_Fluid = (rho * g * Q_ft * h) / 550
    P_Mech = P_Fluid / n_p
    return P_Fluid, P_Mech

def calculate_pulley_tension(P_Mech, RPM, Pulley_rad, T1):
    T2 = ((P_Mech * 550) / ((RPM * 2 * np.pi / 60) * (Pulley_rad / 12))) + T1
    return T2

def calculate_moment_and_torque(T1, T2, L_Pulley_to_Wall, L_Pump_to_Wall):
    M_G = (T1 + T2) * (L_Pulley_to_Wall - L_Pump_to_Wall)
    T_G = (T2 - T1) * (Pulley_rad)
    return M_G, T_G

def solve_forces(T1, T2, W_beam, W_pump, theta, L_Pump_to_Wall, L_support, L_beam, Width_beam, M_G, T_G):
    S1, S2, Wall_X1, Wall_X2, Wall_Z = symbols('S1 S2 Wall_X1 Wall_X2 Wall_Z')
    eq1 = Eq(T1 + T2 - W_pump - S1 * math.sin(theta) - S2 * math.sin(theta) - W_beam + Wall_Z, 0)
    eq2 = Eq(S1 * math.cos(theta) + S2 * math.cos(theta) + (Wall_X1 + Wall_X2), 0)
    eq3 = Eq(M_G + ((T1+T2-W_pump) * L_Pump_to_Wall) - ((S1 * math.sin(theta) + S2 * math.sin(theta)) * L_support) - (W_beam * (L_beam/2)), 0)
    eq4 = Eq(T_G - ((S1*np.sin(theta) + (Wall_Z/2)) * (Width_beam/2))  + ((S2*np.sin(theta) + (Wall_Z/2)) * (Width_beam/2)), 0)
    eq5 = Eq((S1*np.cos(theta) - S2*np.cos(theta) + Wall_X1 - Wall_X2)*0.5*Width_beam, 0)
    soln = solve((eq1, eq2, eq3, eq4, eq5), (S1, S2, Wall_X1, Wall_X2, Wall_Z))

    # print(soln)
    return soln

def calculate_connection_reactions(soln, theta):
    S1_ConnectionX = -soln[symbols('S1')] * math.cos(theta)
    S1_ConnectionZ = -soln[symbols('S1')] * math.sin(theta)
    S2_ConnectionX = -soln[symbols('S2')] * math.cos(theta)
    S2_ConnectionZ = -soln[symbols('S2')] * math.sin(theta)
    return S1_ConnectionX, S1_ConnectionZ, S2_ConnectionX, S2_ConnectionZ


# CONSTRAINTS
Q = 150  # Volume flow rate [GPM]
Q_ft = Q / 7.48 / 60  # Volume flow rate [ft^3/sec]
n_p = 0.65  # Pump Efficiency
h = 60  # Final Height [ft]
rho = 1.94  # slugs/ft3
g = 32.174  # Acceleration Due to Gravity
W_pump = 12.5  # Weight of Pump [lbf]
RPM = 1800  # RPM of Motor
L_Pump_to_Wall = 9.56 + 2.5  # Distance from Pump to Wall [in]
L_Pulley_to_Wall = 18.56  # Distance from Pulley to Wall [in]

# FREE VARIABLES
L_mounting = 6  # Distance between mounting studs [in]
L_support = 5  # Placement [in]
L_beam = 14.56 # Placement [in]
Width_beam = 6  # [in]
Thick_beam = 0.5  # [in]
density_beam = 0.1  # [lb/in^3]

# # FREE VARIABLES
# L_mounting = 6  # Distance between mounting studs [in]
# L_support = 2.00  # Placement [in]
# L_beam = 14.56  # Placement [in]
# Width_beam = 6.00  # [in]
# Thick_beam = 0.38  # [in]

"""
# UNCOMMENT FOR USER INPUT
L_mounting = float(input(print("Enter Distance between mounting studs [in]: ")))
L_support = float(input(print("Enter Distance between support attachment and wall [in]: ")))
L_beam = float(input(print("Enter Length of Beam [in]: ")))
Width_beam = float(input(print("Enter Width of Beam [in]: ")))
Thick_beam = float(input(print("Enter Thickness of Beam [in]: ")))
density_beam = float(input(print("Enter Density of Beam [lb/in^3]: ")))
"""

# Calculated Values
W_beam = density_beam * L_beam * Width_beam * Thick_beam  # Weight of Beam [lbf]
theta = np.arctan(L_mounting / L_support)

P_Fluid, P_Mech = calculate_fluid_mechanical_power(Q_ft, rho, g, h, n_p)
print(f"Fluid Power: {P_Fluid:.6f} hp \nMechanical Power: {P_Mech:.6f} hp")

T1 = 44.9618  # Pulley Tension 1 [lbf]
Pulley_rad = 10 / 2  # [in]
T2 = calculate_pulley_tension(P_Mech, RPM, Pulley_rad, T1)
print(f"Pulley Tension 1: {T1:.6f} lbf")
print(f"Pulley Tension 2: {T2:.6f} lbf")

M_G, T_G = calculate_moment_and_torque(T1, T2, L_Pulley_to_Wall, L_Pump_to_Wall)
print(f"Moment at Point G: {M_G:.6f} lbf-in")
print(f"Torque at Point G: {T_G:.6f} lbf-in")

soln = solve_forces(T1, T2, W_beam, W_pump, theta, L_Pump_to_Wall, L_support, L_beam, Width_beam, M_G, T_G)
print(f"Force S1: {soln[symbols('S1')]:.6f} lbf")
print(f"Force S2: {soln[symbols('S2')]:.6f} lbf")
print(f"Force Wall_X1: {soln[symbols('Wall_X1')]:.6f} lbf")
print(f"Force Wall_X2: {soln[symbols('Wall_X2')]:.6f} lbf")
print(f"Force Wall_Z: {soln[symbols('Wall_Z')]:.6f} lbf")

S1_ConnectionX, S1_ConnectionZ, S2_ConnectionX, S2_ConnectionZ = calculate_connection_reactions(soln, theta)
print(f"Force S1 Wall Connection X: {S1_ConnectionX:.6f} lbf")
print(f"Force S1 Wall Connection Z: {S1_ConnectionZ:.6f} lbf")
print(f"Force S2 Wall Connection X: {S2_ConnectionX:.6f} lbf")
print(f"Force S2 Wall Connection Z: {S2_ConnectionZ:.6f} lbf")


s1 = soln[symbols('S1')]
s2 = soln[symbols('S2')]
wall_X1 = soln[symbols('Wall_X1')]
wall_X2 = soln[symbols('Wall_X2')]
wall_Z = soln[symbols('Wall_Z')]

x_vals_xz = [
    0, 0,
    L_Pulley_to_Wall-L_Pump_to_Wall, L_Pulley_to_Wall-L_Pump_to_Wall,
    L_Pulley_to_Wall-L_beam/2, L_Pulley_to_Wall-L_beam/2,
    L_Pulley_to_Wall - L_support, L_Pulley_to_Wall - L_support,
    L_Pulley_to_Wall, L_Pulley_to_Wall
]

shear_xz = [
    0,
    T1 + T2, T1 + T2,
    (T1 + T2 - W_pump), (T1 + T2 - W_pump),
    (T1 + T2 - W_pump) - W_beam, (T1 + T2 - W_pump) - W_beam,
    (T1 + T2 - W_pump) - W_beam - s1*np.sin(theta)-s2*np.sin(theta), (T1 + T2 - W_pump) - W_beam - s1*np.sin(theta)-s2*np.sin(theta),
    (T1 + T2 - W_pump) - W_beam - s1*np.sin(theta)-s2*np.sin(theta) + wall_Z
]

x_vals_yz = [
    0, 0,
    Pulley_rad-Width_beam/2, Pulley_rad-Width_beam/2,
    Pulley_rad, Pulley_rad,
    Pulley_rad+Width_beam/2, Pulley_rad+Width_beam/2,
    Pulley_rad*2, Pulley_rad*2
]

shear_yz = [
    0,
    T2, T2,
    T2 - s1*np.sin(theta) + wall_Z/2, T2 - s1*np.sin(theta) + wall_Z/2,
    T2 - s1*np.sin(theta) + wall_Z/2 - W_beam - W_pump, T2 - s1*np.sin(theta) + wall_Z/2 - W_beam - W_pump,
    T2 - s1*np.sin(theta) + wall_Z/2 - W_beam - W_pump - s2*np.sin(theta) + wall_Z/2, T2 - s1*np.sin(theta) + wall_Z/2 - W_beam - W_pump - s2*np.sin(theta) + wall_Z/2,
    T2 - s1*np.sin(theta) + wall_Z/2 - W_beam - W_pump - s2*np.sin(theta) + wall_Z/2 + T1
]

unique_x_xz, indices_xz = np.unique(x_vals_xz, return_index=True)
unique_x_yz, indices_yz = np.unique(x_vals_yz, return_index=True)

shear_xz_unique = np.array(shear_xz)[indices_xz]
shear_yz_unique = np.array(shear_yz)[indices_yz]

moment_xz = np.zeros_like(unique_x_xz, dtype=float)
moment_yz = np.zeros_like(unique_x_yz, dtype=float)

moment_xz[0] = shear_xz_unique[0] * unique_x_xz[0]
moment_yz[0] = shear_yz_unique[0] * unique_x_yz[0]

for i in range(1, len(unique_x_xz)):
    dx = unique_x_xz[i] - unique_x_xz[i - 1]
    moment_xz[i] = moment_xz[i - 1] + shear_xz_unique[i] * dx

for i in range(1, len(unique_x_yz)):
    dy = unique_x_yz[i] - unique_x_yz[i - 1]
    moment_yz[i] = moment_yz[i - 1] + shear_yz_unique[i] * dy


fig, axes = plt.subplots(2, 2, figsize=(10, 8))

# XZ Plane Shear Diagram
axes[0, 0].plot(x_vals_xz, shear_xz, marker="o", linestyle="-")
axes[0, 0].axhline(0, color="k", linestyle="--")
axes[0, 0].set_xlabel("X [in]")
axes[0, 0].set_ylabel("V [lbf]")
axes[0, 0].set_title("XZ Plane Shear Diagram")
axes[0, 0].grid(True)

# XZ Plane Moment Diagram
axes[1, 0].plot(unique_x_xz, moment_xz, marker="o", linestyle="-", color="r", label="Moment (M_xz)")
axes[1, 0].axhline(0, color="k", linestyle="--")
axes[1, 0].set_xlabel("X [in]")
axes[1, 0].set_ylabel("M [lbf-in]")
axes[1, 0].set_title("XZ Plane Moment Diagram")
axes[1, 0].grid(True)

# YZ Plane Shear Diagram
axes[0, 1].plot(x_vals_yz, shear_yz, marker="o", linestyle="-")
axes[0, 1].axhline(0, color="k", linestyle="--")
axes[0, 1].set_xlabel("Y [in]")
axes[0, 1].set_ylabel("V [lbf]")
axes[0, 1].set_title("YZ Plane Shear Diagram")
axes[0, 1].grid(True)

# YZ Plane Moment Diagram
axes[1, 1].plot(unique_x_yz, moment_yz, marker="o", linestyle="-", color="r", label="Moment (M_yz)")
axes[1, 1].axhline(0, color="k", linestyle="--")
axes[1, 1].set_xlabel("Y [in]")
axes[1, 1].set_ylabel("M [lbf-in]")
axes[1, 1].set_title("YZ Plane Moment Diagram")
axes[1, 1].grid(True)


# axes[0, 0].set_yticklabels([])
# axes[1, 0].set_yticklabels([])
# axes[0, 1].set_yticklabels([])
# axes[1, 1].set_yticklabels([])

plt.tight_layout()
plt.show()


M_max= np.max(np.abs(moment_xz))
print(f"Peak bending moment in XZ likely at the support connection: (x,y) = ({unique_x_xz[np.argmax(np.abs(moment_xz))]:.3}, {np.max(np.abs(moment_xz)):.6})")
print(f"Peak bending moment in YZ likely in the middle of the beam: (x,y) = ({unique_x_yz[np.argmax(np.abs(moment_yz))]:.3}, {np.max(np.abs(moment_yz)):.6})")

#Platform Stresses:Focused on the Cross Section between Support Connections
Thick_beam= 0.5
Width_beam = 6
L_support = 4.9
#Geometric Constraints
c = Thick_beam/2 #distance to outer fiber

d=0.25 #size of hole
I_EF = (Width_beam/12)*(Thick_beam**3-d**3)
I = (Width_beam*Thick_beam**3)/12
A = (Width_beam-d)*Thick_beam #cross section with pin hole accounted for
k = 0.321

#Material Properties
E = 10000000 #psi for alum 6061
Sy = 40000 #psi for alum 6061

#Stress Concentration Factors
k_bending = 2.1 #Table A-15-2
k_axial = 2.6 #Table A-15-12

#Nominal Stresses for location 1
sigma_bending = (M_max)*(c)/I_EF
sigma_axial = ((S2_ConnectionX+S1_ConnectionX)/A)
torsion = T_G/(k*Width_beam*Thick_beam**2)

#max stresses _ accounting for stress concentrations
sigma_bending_max = k_bending * sigma_bending
sigma_axial_max = k_axial * sigma_axial
torsion_max = torsion
sigma_bending_max = float(sigma_bending_max)
sigma_axial_max = float(sigma_axial_max)
torsion = float(torsion)
sigma1 = 1/np.sqrt(2) * np.sqrt((sigma_bending_max - sigma_axial_max)**2 + 6*(torsion)**2)

print(f"Max Bending Stress for Location 1: {sigma_bending_max:.6f} lbf/in^2")
print(f"Max Axial Stress for Location 1: {sigma_axial_max:.6f} lbf/in^2")
print(f"Max Torsion Stress for Location 1: {torsion_max:.6f} lbf/in^2")
print(f"Max von mises stress at Location 1: {sigma1:.6f} lbf/in2")

#Nominal Stresses for Location 2 at bolt holes  - torsion (neglected at bolt because very small compared to transverse), transverse, axial
V_max = 850
transverse = 575.7  # lbf/in^2
print(f"Max Transverse Shear Stress for Location 2: {transverse:.6f} lbf/in^2")
print(f"Max Axial Stress for Location 2: {sigma_axial_max:.6f} lbf/in^2")

sigma2 = 1/np.sqrt(2) * np.sqrt( sigma_axial_max**2 + 6*(transverse)**2 )
print(f"Max von mises stress at Location 2: {sigma2:.6f} lbf/in^2")

#Max deflection in pump platform
P_G = W_pump-T1-T2
P_EF = S1_ConnectionZ+S2_ConnectionZ #CHECK SIGN

delta_G = ((P_G*L_Pump_to_Wall**2)*(3*L_beam-L_Pump_to_Wall))/(6*E*I) + (M_G*L_beam**2)/(2*E*I) #deflection at the tip of platform due to forces at F
delta_supports = ((P_EF*L_support**2)*(3*L_beam-L_support))/(6*E*I) #deflection at the tip of platform due to forces from the supports
delta_max = delta_G + delta_supports #assuming the weight of the platform has negligible effects on the deflection
print(f"Max Deflection in Pump Platform: {delta_max:.6f} in") #negative deflection means the platform deflects upwards

#Max shear and principal stress calcs
C = (sigma_bending_max+sigma_axial_max)/2
R= np.sqrt(((sigma_bending_max+sigma_axial_max)/2)**2 + torsion**2)
shear_max = R
principal_max = C + R
print(f"Max Shear for Location 1: {shear_max:.6f} lbf/in^2")
print(f"Max Principal Stress for Location 1: {principal_max:.6f} lbf/in^2")

#Factor of Safety for Platform
FoS = Sy/(sigma1) #ductile material so ductile criterion
print(f"Factor of Safety for Platform: {FoS:.6f}")


#Support Arm Stresses

#Geometric Constraints
#Mounting hole 0 - 3in, 1 - 9in, 2 - 15in, 3 - 21in
mounting_hole = 2
Support_arm_length = math.sqrt((mounting_hole*L_mounting)**2+L_support**2)
Support_arm_thick = 0.25
Support_arm_width = 0.75
Support_arm_effective_l_pinned = 1.0  #pinned-pinned
Support_arm_effective_l_fixed = 0.5 #fixed-fixed
Support_arm_area = (Support_arm_width - d) * Support_arm_thick  #accounting for pin hole
Support_arm_I_xx = (Support_arm_width**3 * Support_arm_thick) / 12
Support_arm_I_yy = (Support_arm_width * Support_arm_thick**3) / 12

#Force in Support Arm B-F
S1 = soln[symbols('S1')]

#Nominal Stress
Support_arm_sigma_0 = S1 / Support_arm_area

#Stress Concentration Factors
Support_arm_d_w_ratio = d / Support_arm_width
Support_arm_L_w_ratio = (Support_arm_length - Support_arm_width / 2) / Support_arm_width
print(f"Support arm length: {Support_arm_length:.6f}")
print(f"d/w ratio: {Support_arm_d_w_ratio:.6f}")
print(f"L/w ratio: {Support_arm_L_w_ratio:.6f}")
Support_arm_K_t = 3.5 #Table A-15-12

#Maximum Stresses accounting for Stress Concentration
Support_arm_sigma_x = Support_arm_K_t * Support_arm_sigma_0

#Stress Concentration
sig_xx_BF = Support_arm_sigma_x
sig_yy_BF = 0
tau_xy_BF = 0

#Mohr's Circle
CMC_BF = (sig_xx_BF + sig_yy_BF)/2
RMC_BF = math.sqrt(0.5*(sig_xx_BF-sig_yy_BF)**2 + 0.5*tau_xy_BF**2)
sig1_BF = CMC_BF + RMC_BF
sig2_BF = CMC_BF - RMC_BF
tauMax_BF = RMC_BF

#Ductile Failure (von Mises Distortion Energy)
sig_von_BF = math.sqrt(sig1_BF**2 + sig2_BF**2 - sig1_BF*sig2_BF)
#Factor of safety
n_BF = Sy/sig_von_BF

# Buckling
C_arm2 = 1
b_BF = min(Support_arm_thick, Support_arm_width)
h_BF = max(Support_arm_thick, Support_arm_width)
I_BF = (b_BF*h_BF**3)/12
A_BF = Support_arm_thick*Support_arm_width
k_BF = math.sqrt(I_BF/A_BF)
l_BF = Support_arm_length
lk_BF_1 = ((2*(math.pi**2)*C_arm2*E)/Sy)**0.5

# Euler buckling vs. JB Johnson buckling
if(l_BF/k_BF >= lk_BF_1):
  Pcr_BF = b_BF*h_BF*C_arm2*(math.pi**2)*E/(l_BF/(h_BF/math.sqrt(12)))**2
else:
  Pcr_BF = Sy - (1/(C_arm2*E))*(Sy/(2*math.pi*(h_BF/math.sqrt(12))))**2

#Factors of Safety
#Critical buckling is larger than sigma_1
Support_arm_n_buckling = Pcr_BF / S1

#Maximum Deflection
Support_arm_delta_max = (S1 * Support_arm_length) / (Support_arm_area * E)

print(f"Nominal Stress: {Support_arm_sigma_0:.6f} lbf/in^2")
print(f"Max Stress (with K_t): {Support_arm_sigma_x:.6f} lbf/in^2")
print(f"Max Shear Stress: {tauMax_BF:.6f} lbf/in^2")
print(f"Principal Stresses: Sigma_1 = {sig1_BF:.6f}, Sigma_2 = {sig2_BF:.6f}")
print(f"Critical Buckling Load: {Pcr_BF:.6f} lbf")
print(f"Factor of Safety Buckling: {Support_arm_n_buckling:.6f}")
print(f"Factor of Safety (Von Mises): {n_BF:.6f}")
print(f"Max Deflection: {Support_arm_delta_max:.6f} in")

#pin in support
#dimensions
r = 0.120
A =  np.pi * r**2

#material properties
Sy = 51000 #yield strength for 6061 alum in PSI
#shear force
V_support = soln[symbols('S1')]

#shear stress calcs
tau_support = V_support/A
print(f"Max Shear Stress: {tau_support:.6f} lbf/in^2")

#pin at wall
V_wall = np.sqrt((169.237966)**2+(875.941312)**2)
print(V_wall)
tau_wall = V_wall/A
print(f"Max Shear Stress: {tau_wall:.6f} lbf/in^2")


Fos_SupportPin = 0.577*Sy/(tau_support)
print(f"Factor of Safety for Pin in Support: {Fos_SupportPin:.6f}")
Fos_WallPin = 0.577*Sy/(tau_wall)
print(f"Factor of Safety for Pin in Wall: {Fos_WallPin:.6f}")

# Define your data
data = {
    "Trial": [1, 2, 3, 4, 5, 6, 7, 8],
    "# Elements": [3312, 4213, 5211, 7794, 10797, 11300, 13499, 19771],
    "Platform (psi)": [10330, 5879, 6560, 8405, 6947, 6925, 6325, 6676],
    "Arm 1 (psi)": [7314, 3798, 4378, 6872, 5098, 5293, 4516, 4943],
    "Arm 2 (psi)": [6557, 3392, 3653, 5863, 4095, 4086, 3708, 4076],
    "Deflection (mm)": [4.666, 2.802, 3.01, 3.337, 3.232, 3.254, 3.059, 3.14]
}

df = pd.DataFrame(data)

# Parts to plot
parts = {
    "Platform (psi)": "navy",
    "Arm 1 (psi)": "navy",
    "Arm 2 (psi)": "navy"
}

# Plot each part
for part, color in parts.items():
    x = np.array(df["# Elements"])
    y = np.array(df[part])
    y_lower = y * 0.95
    y_upper = y * 1.05

    if part == "Platform (psi)":
        fig, ax1 = plt.subplots(figsize=(8, 5))

        # Stress line and min/max bounds
        ax1.plot(x, y, marker='o', linestyle='-', color=color, label=part)
        for xi, yi_low, yi_high in zip(x, y_lower, y_upper):
            ax1.plot([xi, xi], [yi_low, yi_high], color='gray', linewidth=2)  # Vertical range line

        ax1.set_xlabel("Number of Elements")
        ax1.set_ylabel("Stress (psi)", color=color)
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.set_title("Mesh Convergence - Platform")
        ax1.grid(True)

        # Deflection on secondary axis
        ax2 = ax1.twinx()
        ax2.plot(x, df["Deflection (mm)"], marker='D', linestyle='--', color='red', label='Deflection (mm)')
        ax2.set_ylabel("Deflection (mm)", color='red')
        ax2.tick_params(axis='y', labelcolor='red')

        # Combine legends
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc='center right')

        plt.tight_layout()
        plt.show()

    else:
        plt.figure(figsize=(8, 5))
        plt.plot(x, y, marker='o', linestyle='-', color=color, label=part)
        for xi, yi_low, yi_high in zip(x, y_lower, y_upper):
            plt.plot([xi, xi], [yi_low, yi_high], color='gray', linewidth=2)

        plt.xlabel("Number of Elements")
        plt.ylabel("Stress (psi)")
        plt.title(f"Mesh Convergence - {part.split()[0]} {part.split()[1]}")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()
