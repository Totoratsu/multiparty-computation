# MPC from Somewhat Homomorphic Encryption (DPSZ Protocol)

Este repositorio contiene una implementación en **Python** y **SageMath** del protocolo de **Computación Multipartita Segura (MPC)** basado en Cifrado Algo Homomórfico (SHE), siguiendo la arquitectura descrita en el paper *Multiparty Computation from Somewhat Homomorphic Encryption* (Damgård, Pastro, Smart, Zakarias - 2011).

El proyecto implementa la **Fase de Preprocesamiento (Offline Phase)**, la cual genera material criptográfico correlacionado (Triples de Beaver, Pares aleatorios) autenticados mediante MACs, utilizando una variante del esquema de cifrado **BGV** con operaciones SIMD.

## Características Principales

  * **Esquema de Cifrado BGV:** Implementación sobre anillos de polinomios $\mathcal{R}_q$ y $\mathcal{R}_p$ con soporte para operaciones SIMD (Single Instruction, Multiple Data).
  * **Zero-Knowledge Proofs (ZKPoPK):** Protocolo Sigma no interactivo (transformación Fiat-Shamir) para demostrar conocimiento del texto plano y aleatoriedad sin revelar secretos.
  * **Sub-protocolos de MPC:**
      * `Reshare`: Conversión de un cifrado público en participaciones (shares) aditivas privadas.
      * `PBracket` ($\Pi_{Bracket}$): Generación de valores compartidos con autenticación MAC.
      * `PAngle` ($\Pi_{Angle}$): Generación de claves globales y MACs sobre valores cifrados.
  * **Fase de Preprocesamiento Completa ($\Pi_{PREP}$):**
      * `Initialize`: Generación distribuida de la clave MAC global $\alpha$.
      * `Pair`: Generación de pares aleatorios $(r, [r], \langle r \rangle)$.
      * `Triple`: Generación de Triples de Beaver $(a, b, c)$ autenticados.

## Estructura del Proyecto

```text
.
├── preprocessing/          # Módulos del núcleo criptográfico
│   ├── encryption.py       # Funciones BGV (encrypt, decrypt, homomorphic ops)
│   ├── zk_protocols.py     # Pruebas de conocimiento cero (ZKPoPK)
│   └── functions.py        # Utilidades (muestreo gaussiano, CRT, encoding)
├── protocols/              # Lógica de alto nivel del protocolo MPC
│   └── preprocessing.py    # Clase PreprocessingProtocol (Initialize, Pair, Triple)
├── tests/                  # Scripts de prueba y verificación
└── README.md
```

## Desarrollo

Este proyecto depende fuertemente de **SageMath** para la aritmética de campos finitos y anillos de polinomios. Se recomienda utilizar Conda para gestionar el entorno.

### Requisitos Previos

  * [Conda](https://docs.conda.io/en/latest/) (Miniconda o Anaconda)

### Configuración del Entorno

Para instalar las dependencias y configurar el entorno de desarrollo, ejecuta:

```bash
# Crear el entorno con Python 3.12 y SageMath desde conda-forge
conda create -n sage-env python=3.12 sage -c conda-forge

# Activar el entorno
conda activate sage-env
```

## Uso

### Ejecutar Pruebas

El proyecto incluye un script de prueba integral que verifica la corrección matemática de los triples generados ($c = a \cdot b$) y la integridad del descifrado distribuido.

```python
from protocols.preprocessing import PreprocessingProtocol
# ... (configuración de contexto) ...

# Instanciar el protocolo
proto = PreprocessingProtocol(context)

# 1. Inicializar claves MAC
proto.run_initialize()

# 2. Generar un Triple de Beaver
a_shares, b_shares, c_shares, angle_a, angle_b, angle_c = proto.run_triple()

# Los resultados (a_shares, etc.) son vectores aditivos distribuidos entre los jugadores.
```

### Parámetros de Seguridad

La implementación permite configurar los parámetros del sistema BGV:

  * `N`: Grado del polinomio ciclotómico (dimensión del retículo).
  * `q`: Módulo del texto cifrado (debe ser suficientemente grande para soportar el ruido).
  * `p`: Módulo del texto plano.
  * `s`: Número de slots SIMD.
  * `sigma`: Desviación estándar para el ruido gaussiano discreto.

## Referencias

  * **[DPSZ11]** Damgård, I., Pastro, V., Smart, N., & Zakarias, S. (2011). *Multiparty Computation from Somewhat Homomorphic Encryption*. CRYPTO 2012. [Enlace al paper](https://eprint.iacr.org/2011/535.pdf).
  * **[BGV11]** Brakerski, Z., Gentry, C., & Vaikuntanathan, V. (2011). *Fully Homomorphic Encryption without Bootstrapping*.