#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <sstream>
#include <vector>

using namespace std;

// ================= PARÂMETROS DE REFERÊNCIA =================
const double T_REF = 392.41;       // Temperatura da simulação base (K)
const double DP_REF = 4550.0;      // Queda de pressão base (Pa) -> 4.55 kPa
const double W1 = 0.5;             // Peso para a Temperatura
const double W2 = 0.5;             // Peso para a Queda de Pressão
// ============================================================

// Função para extrair o último valor numérico de um arquivo .dat do OpenFOAM
double lerUltimoValorDat(const string& caminhoArquivo) {
    ifstream file(caminhoArquivo);
    if (!file.is_open()) {
        cerr << "Erro: Nao foi possivel abrir o arquivo " << caminhoArquivo << endl;
        return 1e9; // Retorna um valor alto para penalizar a falha
    }

    string linha, ultimaLinhaValida;
    while (getline(file, linha)) {
        // Ignora linhas de cabeçalho (que começam com #) e linhas vazias
        if (!linha.empty() && linha[0] != '#') {
            ultimaLinhaValida = linha;
        }
    }
    file.close();

    // Extrai o valor da segunda coluna (ignora a coluna do tempo/iteração)
    stringstream ss(ultimaLinhaValida);
    double tempo, valor = 1e9;
    ss >> tempo >> valor;
    
    return valor;
}

// Função para calcular a Função Objetivo multiobjetivo
double calcularFuncaoObjetivo(double vel_mag) {
    cout << "\n---------------------------------------------------" << endl;
    cout << "[C++] Testando velocidade de injeção: " << vel_mag << " m/s" << endl;

    // 1. Definição do vetor velocidade (Injeção fixa no eixo Z)
    double Ux = 0.0;
    double Uy = 0.0;
    double Uz = vel_mag; // Toda a magnitude direcionada no eixo do canal

    // 2. Manipulação do arquivo 0/fluid/U
    ofstream fileU("0/fluid/U");
    fileU << "/*--------------------------------*- C++ -*----------------------------------*\\\n";
    fileU << "| =========                 |                                                 |\n";
    fileU << "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n";
    fileU << "|  \\\\    /   O peration     | Version:  v2412                                 |\n";
    fileU << "|   \\\\  /    A nd           | Website:  www.openfoam.com                      |\n";
    fileU << "|    \\\\/     M anipulation  |                                                 |\n";
    fileU << "\\*---------------------------------------------------------------------------*/\n";
    fileU << "FoamFile { version 2.0; format ascii; class volVectorField; object U; }\n";
    fileU << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";
    fileU << "dimensions      [0 1 -1 0 0 0 0];\n";
    fileU << "internalField   uniform (0 0 0.1);\n\n";
    fileU << "boundaryField\n{\n";
    fileU << "    inlet\n    {\n";
    fileU << "        type            fixedValue;\n";
    fileU << "        value           uniform (" << Ux << " " << Uy << " " << Uz << ");\n";
    fileU << "    }\n";
    fileU << "    outlet { type inletOutlet; inletValue uniform (0 0 0); value uniform (0 0 0); }\n";
    fileU << "    bobinas { type noSlip; }\n";
    fileU << "    walls { type noSlip; }\n";
    fileU << "}\n";
    fileU << "// ************************************************************************* //\n";
    fileU.close();

    // 3. Execução automatizada (Paralelo, Multirregião e Limpeza Segura)
    cout << "[C++] Limpando pastas antigas..." << endl;
    system("rm -rf processor*");             
    system("foamListTimes -rm 2>/dev/null"); 
    system("rm -rf postProcessing");         

    cout << "[C++] Decompondo o dominio (aplicando nova velocidade nos nucleos)..." << endl;
    system("LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6 decomposePar -allRegions -force > log.decompose");

    cout << "[C++] Rodando chtMultiRegionSimpleFoam em PARALELO..." << endl;
    int status = system("LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6 mpirun -np 20 chtMultiRegionSimpleFoam -parallel > log.foam");

    // Trava de segurança
    if (status != 0) {
        cout << "[!] ERRO FATAL: A simulacao falhou ou foi abortada manualmente!" << endl;
        return 1e12; 
    }

    // 4. Pós-Processamento
    string pathTMedia = "postProcessing/fluid/tempInterfaceBobina_Media/0/surfaceFieldValue.dat";
    string pathPressaoIn = "postProcessing/fluid/quedaPressaoInlet/0/surfaceFieldValue.dat";
    string pathPressaoOut = "postProcessing/fluid/quedaPressaoOutlet/0/surfaceFieldValue.dat";

    double T_media = lerUltimoValorDat(pathTMedia);
    double P_in = lerUltimoValorDat(pathPressaoIn);
    double P_out = lerUltimoValorDat(pathPressaoOut);
    
    double delta_P = abs(P_in - P_out);

    // 5. Cálculo da Função Objetivo
    double f_obj = W1 * (T_media / T_REF) + W2 * (delta_P / DP_REF);
    
    cout << "[C++] Resultados para " << vel_mag << " m/s:" << endl;
    cout << "      T_media = " << T_media << " K" << endl;
    cout << "      Delta P = " << delta_P << " Pa" << endl;
    cout << "      F(obj)  = " << f_obj << endl;

    return f_obj;
}

int main() {
    cout << "Iniciando Otimizacao via Razao Aurea (Velocidade)..." << endl;

    // Limites de busca estreitados para acelerar a convergencia
    double a = 0.3; 
    double b = 1.3; 
    
    const double phi = (sqrt(5.0) - 1.0) / 2.0; 
    
    // Tolerancia aumentada para reduzir os passos
    double tolerancia = 0.1; 

    double c = b - phi * (b - a);
    double d = a + phi * (b - a);

    double fc = calcularFuncaoObjetivo(c);
    double fd = calcularFuncaoObjetivo(d);

    int iteracao = 1;

    while ((b - a) > tolerancia) {
        cout << "\n--- Iteracao " << iteracao << " ---" << endl;
        cout << "Intervalo atual: [" << a << " , " << b << "] m/s" << endl;

        if (fc < fd) {
            b = d;
            d = c;
            fd = fc;
            c = b - phi * (b - a);
            fc = calcularFuncaoObjetivo(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + phi * (b - a);
            fd = calcularFuncaoObjetivo(d);
        }
        iteracao++;
    }

    double velocidade_otima = (a + b) / 2.0;
    cout << "\n=================================================" << endl;
    cout << "OTIMIZACAO CONCLUIDA!" << endl;
    cout << "Velocidade ideal estimada: " << velocidade_otima << " m/s" << endl;
    cout << "=================================================" << endl;

    return 0;
}
