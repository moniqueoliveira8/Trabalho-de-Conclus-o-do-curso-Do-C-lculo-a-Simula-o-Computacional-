# Otimização Fluidodinâmica de Canais de Resfriamento pelo Método da Razão Áurea e Interpolação 3D de Perdas Eletromagnéticas via Newton, Splines e RBF

Projeto desenvolvido como parte do curso **Do Cálculo à Simulação Computacional** (L2C — Soluções em Computação Científica), sob orientação do **Prof. Dr. Rafael Gabler Gontijo**.

**Autora:** Monique Caroline Silva Oliveira

---

## Sobre o projeto

Este trabalho aborda o gerenciamento térmico de um motor elétrico de fluxo axial com arquitetura **YASA** (*Yokeless and Segmented Armature*, 12 bobinas / 10 polos), resfriado por um circuito de óleo dielétrico. O motor foi modelado inteiramente de forma autoral em **FreeCAD**, servindo de base geométrica para todas as simulações.

Antes de qualquer otimização, o modelo foi validado com base na formulação de **Wanjiku et al. (2021)**, que unifica as análises térmica e elétrica para essa topologia de máquina. A partir dessa base validada, o trabalho resolve dois problemas de engenharia complementares:

1. **Otimização fluidodinâmica do canal de resfriamento**, via **Método da Razão Áurea**, buscando a velocidade de injeção do óleo que minimiza a temperatura das bobinas sem penalizar excessivamente a perda de carga hidráulica.
2. **Interpolação 3D das perdas eletromagnéticas**, transferindo o campo de perdas calculado no **ANSYS Maxwell** (malha de Elementos Finitos) para a malha de **Volumes Finitos** do **OpenFOAM**, usando a expansão espacial das Diferenças Divididas de Newton / *splines* (via Tetraedrização de Delaunay) combinada com Funções de Base Radial (RBF) nas regiões de fronteira.

## Problema 1 — Otimização via Razão Áurea

- Duas geometrias candidatas de canal (*single channel 1* e *single channel 2*) foram avaliadas via CFD conjugado (*solver* CHT, modelo de turbulência $k$-$\omega$ SST) em malhas de ~7 milhões de elementos.
- Como o canal de admissão em formato de "T" mitiga a influência do ângulo de injeção, a variável de otimização foi redefinida para a **magnitude da velocidade de injeção**.
- Um algoritmo em **C++** automatiza todo o ciclo: define o intervalo de busca segundo a proporção áurea ($\phi \approx 0{,}618$), edita o dicionário `0/fluid/U` do OpenFOAM, decompõe e executa o *solver* `chtMultiRegionSimpleFoam` em paralelo (MPI), lê os arquivos `.dat` de pós-processamento e calcula uma função objetivo ponderada entre temperatura média das bobinas e perda de carga ($\Delta P$).
- A cada iteração, o subintervalo de pior desempenho é descartado, convergindo para a velocidade ótima de injeção.
- Resultado: velocidade ótima convergida de **1,25492 m/s**, com $F_{obj} = 0{,}883995$.

## Problema 2 — Interpolação Eletromagnética-Térmica

- As perdas eletromagnéticas (densidade de fluxo magnético $B_{pico}$) são exportadas do ANSYS Maxwell (FEM) e convertidas em perdas Joule contínuas (DC) e por correntes parasitas (AC), consolidadas no termo fonte $Q_{em}$ [W/m³].
- Como as malhas FEM (ANSYS) e FVM (OpenFOAM) não coincidem espacialmente, a transferência é feita por um algoritmo em **Python**:
  - **Prova de conceito 1D**: validação do formalismo de Diferenças Divididas de Newton e *splines* em malhas lineares, antes da expansão 3D.
  - **Tetraedrização de Delaunay + `LinearNDInterpolator`**: preenchimento do volume interno das bobinas (~2 milhões de elementos na malha térmica).
  - **RBF (`RBFInterpolator`)**: aplicada nas células de fronteira/casca externa, restrita aos 30 vizinhos mais próximos, garantindo continuidade $C^1/C^2$ e evitando valores nulos (NaN) ou saltos não físicos que fariam o *solver* térmico divergir.
- O campo final é exportado como `volScalarField` (`constant/Qem`), pronto para a simulação térmica conjugada no OpenFOAM.

## Ferramentas e tecnologias

- **FreeCAD** — modelagem geométrica autoral do motor
- **ANSYS Maxwell** — simulação eletromagnética (FEM)
- **OpenFOAM** — CFD e transferência de calor conjugada (FVM, `chtMultiRegionSimpleFoam`)
- **C++** — automação da otimização via Razão Áurea
- **Python** (`scipy.interpolate`: `LinearNDInterpolator`, `RBFInterpolator`) — interpolação 3D dos dados eletromagnéticos
- **ParaView** — pós-processamento e visualização

## Estrutura do repositório

```
.
├── simulacao_openfoam/       # Casos e configurações do OpenFOAM
├── otimizador_cpp/           # Algoritmo de automação da Razão Áurea
├── interpolacao_python/      # Scripts de mapeamento Newton/Delaunay + RBF
├── logs/                     # Históricos de execução dos ensaios de otimização
└── docs/                     # Relatório completo e figuras
```

## Status

Etapa eletromagnética-fluidodinâmica concluída. Próximos passos: rodar a otimização com o número completo de iterações necessário para a convergência térmica em regime permanente (~4.000 iterações) e explorar variações geométricas do canal de admissão.

## Referências principais

- WANJIKU, J. et al. *Electromagnetic and direct-cooling analysis of a traction motor*. IEEE ECCE, 2021.
- TARAN, N. et al. *An overview of methods and a new 3D FEA and analytical hybrid technique for calculating AC winding losses in PM machines*.
- MENTER, F. R. *Two-equation eddy-viscosity turbulence models for engineering applications*. AIAA Journal, 1994.
- GONTIJO, R. G. *Manual de Cálculo Numérico Aplicado*. Eduzz.
