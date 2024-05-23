import numpy as np
import cantera as ct
import argparse

parser = argparse.ArgumentParser(
    description='Generate Lewis numbers for a given mechanism using a 1D premixed flame simulated in Cantera'
)

parser.add_argument('mechanism', type=str, help='Path to the mechanism file')
parser.add_argument('-plot',     action='store_true', help='Flag to plot Lewis numbers across the domain')
parser.add_argument('--output',  type=str, default="LewisNumbersDict", help='Path to the output file')

parser.add_argument('-Lew',         type=float, default=0.5, help='Weight of the Lewis number of the burnt part (0 to 1). 0 - fresh gas, 1 - burnt gas, 0.5 - average. Default: 0.5')
parser.add_argument('-p',           type=float, default=101325, help='Pressure in Pa')
parser.add_argument('-T',           type=float, default=300.0, help='Unburned gas temperature in K')
parser.add_argument('-w',           type=float, default=0.03, help='Width of the domain in m') # DUPLICATION
parser.add_argument('-ER',          type=float, default=1.0, help='Equivalence ratio')
parser.add_argument('-fuel',        type=str,   default="CH4", help='Fuel')
parser.add_argument('-oxidizer',    type=str,   default="O2:0.21,N2:0.79", help='Oxidizer')
parser.add_argument('-transport',   type=str,   default='Multi', help='Transport model')
args = parser.parse_args()

path_img  = 'Le_1d_premixed.png'

print("Selected parameters:")
print(f"Pressure = {args.p}")
print(f"Unburned gas temperature = {args.T}")
print(f"Width of the domain = {args.w}")
print(f"Equivalence ratio = {args.ER}")
print(f"Fuel = {args.fuel}")
print(f"Oxidizer = {args.oxidizer}")
print(f"Transport model = {args.transport}")
print()

gas = ct.Solution(args.mechanism)
gas.TP = args.T, args.p
gas.set_equivalence_ratio(args.ER, args.fuel, args.oxidizer)
f = ct.FreeFlame(gas, width=args.w)
f.transport_model = args.transport
f.energy_enabled = True

f.solve(0)


thermal_diff = f.thermal_conductivity/f.cp_mass/f.density_mass

if args.plot:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=[12,8])
fo = open(args.output, "w")
print_entry = "{:8} {:5.2f}          {:1.5f}         {:7.2f}"
write_entry = "{:8} {};\n"
print('Specie     Le   Delta b/w fresh and burnt  Rel deviation [%]')
for i, specie in enumerate(gas.species_names):
    Lewis_flame = thermal_diff / f.mix_diff_coeffs[i]
    Le = (1-args.Lew)*Lewis_flame[0] + args.Lew*Lewis_flame[-1]
    
    delta_Le = np.abs(Lewis_flame[0]-Lewis_flame[-1])
    print(print_entry.format(specie, Le, delta_Le, delta_Le/Le*100))
    fo.write(write_entry.format(specie, Le))

    if args.plot:
        ax.plot(f.grid, Lewis_flame, label=specie)

fo.close()
print('Lewis numbers were written to ' + args.output)

if args.plot:
    ax.legend(ncol=2)
    ax.set_xlabel('x, m')
    ax.set_ylabel('Le')
    fig.savefig(path_img)
    print('Validation image of 1d premixed flame was written to ' + path_img)
