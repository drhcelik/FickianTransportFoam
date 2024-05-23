import math
import numpy as np
import cantera as ct

import argparse
parser = argparse.ArgumentParser(
    description='Get binary diffusion coefficients for a given mechanism from Cantera and write them to a file'
)
parser.add_argument('mechanism', type=str, help='Path to the mechanism file')
parser.add_argument('-o', '--output', type=str, default="binaryDiffusionCoefficientsPolynomialDict", help='Path to the output file')

args = parser.parse_args()

gas = ct.Solution(args.mechanism)

N_SP = len(gas.species_names)
entry = """{} {{
    type binaryDiffusionCoefficientsPolynomial;
    polynomialCoeffs ({});
}}
"""

f = open(args.output, "w")

for i in range(N_SP):
    for j in range(N_SP):
        if (j>i):
            continue
        pair = gas.species_names[i] + "-" + gas.species_names[j]
        coeff = str(gas.get_binary_diff_coeffs_polynomial(i,j))[1:-1]
        f.write(entry.format(pair, coeff))

f.close()
print('Binary diffusion coefficients were written to binaryDiffusionCoefficientsPolynomialDict')
