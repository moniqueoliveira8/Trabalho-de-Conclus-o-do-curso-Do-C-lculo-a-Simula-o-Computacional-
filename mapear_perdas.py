import numpy as np
import pandas as pd
import re
import os
from scipy.interpolate import LinearNDInterpolator, RBFInterpolator

print("=========================================================================")
print(" MAPEAMENTO MULTIFÍSICO: PERDAS AC/DC E INTERPOLAÇÃO 3D (NEWTON/RBF)")
print("=========================================================================\n")

print("1. Extração da Dispersão de Dados (Nó a Nó)...")
nome_arquivo = 'cindy.csv' if os.path.exists('cindy.csv') else 'cindy.cvs'

# O conjunto de dados originais equivale aos n+1 pontos a serem ajustados
df_em = pd.read_csv(nome_arquivo)
pontos_em = df_em[['Points:0', 'Points:1', 'Points:2']].values
valores_B = df_em['Mag_B'].values 

# ---------------------------------------------------------------------
# PARÂMETROS DA BOBINA
# ---------------------------------------------------------------------

print("-> Calculando Perdas Eletromagnéticas (AC + DC)...")
# Parâmetros geométricos e materiais da bobina
w = 0.0025       # Largura do condutor individual (m) [2.5 mm]
h = 0.0050       # Altura do condutor individual (m) [5.0 mm]
omega = 376.99   # Frequência angular elétrica (rad/s) = 2 * pi * 60 Hz
sigma = 5.8e7    # Condutividade elétrica do cobre a 20°C (S/m)

# Parâmetros macroscópicos para o cálculo DC
l = 0.15         # Comprimento ativo do condutor (m) [150 mm]
Nc = 2           # Número de condutores paralelos
Nt = 20          # Número de espiras por bobina
Ns = 12          # Número de bobinas no estator
I_rms = 15.0     # Corrente RMS do condutor (A)

# ---------------------------------------------------------------------
# CÁLCULO DA PERDA AC (Eddy Currents Volumétricas ponto a ponto)
# ---------------------------------------------------------------------
# Unidade resultante mantida estritamente em W/m³
termo_fracao_ac = ((omega**2) * sigma) / 24.0
termo_geom = (w**2) + (h**2)

valores_perda_ac = termo_fracao_ac * ((valores_B**2)) * termo_geom

# ---------------------------------------------------------------------
# CÁLCULO DA PERDA DC (Constante Volumétrica)
# ---------------------------------------------------------------------
area_condutor = w * h
comprimento_total = l * Nc * Nt * Ns

# 1. Resistência DC da bobina (R = L / (sigma * A))
R_dc = comprimento_total / (sigma * area_condutor)

# 2. Perda Total DC em Watts (P = R * I^2)
P_dc_total = R_dc * (I_rms**2)

# 3. Volume total da bobina
volume_bobina = comprimento_total * area_condutor

# 4. Perda DC Volumétrica Constante em W/m³
P_dc_vol = P_dc_total / volume_bobina
print(f"   [INFO] Perda DC Constante calculada: {P_dc_vol:.4f} W/m³")

# ---------------------------------------------------------------------
# SOMA TOTAL (Qem final)
# ---------------------------------------------------------------------
# Perda contínua (DC) + Gradiente (AC)
valores_perda_total = valores_perda_ac + P_dc_vol
# =====================================================================

print("\n2. Leitura do Domínio Discretizado (OpenFOAM)...")
def extrair_pontos_foam(caminho_arquivo):
    pontos = []
    with open(caminho_arquivo, 'r') as f:
        conteudo = f.read()
        bloco = re.search(r'internalField\s+nonuniform\s+List<vector>\s+\d+\s*\((.*?)\)\s*;', conteudo, re.DOTALL)
        if bloco:
            linhas = bloco.group(1).strip().split('\n')
            for linha in linhas:
                linha_limpa = linha.replace('(','').replace(')','').strip()
                if linha_limpa:
                    pontos.append([float(x) for x in linha_limpa.split()])
    return np.array(pontos)

pontos_foam = extrair_pontos_foam('constant/C')
n_cells = len(pontos_foam)
print(f"   -> Domínio lido: {n_cells} posições-alvo para interpolação.")

print("\n3. Execução dos Procedimentos de Interpolação Numérica...")
# Fundamento 1: Extensão 3D do conceito de Splines Lineares
print("   -> Aplicando Spline Linear N-Dimensional (Tetraedrização de Delaunay)...")
interpolador_linear = LinearNDInterpolator(pontos_em, valores_perda_total)
valores_interpolados = interpolador_linear(pontos_foam)

# Identificação de pontos fora da malha base (necessidade de extrapolação)
valores_nan = np.isnan(valores_interpolados)

if valores_nan.any():
    num_nan = np.sum(valores_nan)
    # Fundamento 2: Continuidade de Derivadas (Análogo a Splines Superiores)
    print(f"   -> Detectados {num_nan} pontos na fronteira geométrica exterior.")
    print("   -> Aplicando Funções de Base Radial (RBF) para impor continuidade do gradiente...")
    
    interpolador_rbf = RBFInterpolator(pontos_em, valores_perda_total, neighbors=30)
    valores_interpolados[valores_nan] = interpolador_rbf(pontos_foam[valores_nan])

print("\n4. Exportando VolScalarField Numérico (constant/Qem)...")
header = f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           |
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "constant";
    object      Qem;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [1 -1 -3 0 0 0 0]; // Unidade de W/m³

internalField   nonuniform List<scalar> {n_cells}
(
"""

footer = """
);

boundaryField
{
    ".*"
    {
        type            zeroGradient;
    }
}

// ************************************************************************* //
"""

with open('constant/Qem', 'w') as f:
    f.write(header)
    for v in valores_interpolados:
        f.write(f"{v}\n")
    f.write(footer)

print("\n[SUCESSO] Malha processada e equacionada perfeitamente. Arquivo gerado.")
