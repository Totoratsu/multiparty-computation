import gmpy2
import numpy as np
from sage.all import Zmod, PolynomialRing, ZZ, GF

from .A_q_space import calculate_m, generate_Aq
from .monte_carlo_C_m import calculate_Cm_parallel, unique_prime_factors
from .fixed_point_params import fixed_point_qr


################################################################
#                    ParamGen(1^k, M)                          #
################################################################

def ParamGen(p, k, s, n, sec):
    # General parameters
    tau = p/2
    delta = 1.0052

    # Generate m, phi_m
    min_m = 14620
    max_m = 50000  # rango de búsqueda de índices ciclotómicos m
    m, phi_m, coeffs = calculate_m(p, k, s, min_m, max_m)
    N = phi_m

    c_sec = 9 * np.power(float(N), 2) * np.power(float(sec), 4) * np.power(2.0, sec + 8)
    
    # Generate G_m
    factors = unique_prime_factors(m)
    if len(factors) == 1:
        C_m = 4/np.pi
        print("m es primo, asi que C_m =", C_m)
    elif np.prod(factors) <= 400:
        C_m = 8.6
        print("los factores de m son pequeños, entonces C_m =", C_m)
    else:
        # Ejecutar en paralelo (aumenta trials para mayor precisión)
        C_m = float(calculate_Cm_parallel(m, total_trials=1000))

        print(f"\nResultado Monte Carlo Paralelo:")
        print(f"C_m estimado: {C_m:.2f}")

        # Margen de seguridad (x1.5 o x2 es común en la práctica)
        margin = 1
        cm_seguro = C_m * margin

        print(f"Valor recomendado (con margen x{margin}): {cm_seguro:.2f}")

    # Generate q, r
    r_seed = 3.2 

    # Ejecutar el método numérico
    q, r = fixed_point_qr(r_seed, N, C_m, p, tau, c_sec, n, sec, delta)

    print("\n=============================================")
    print("      RESULTADOS DE PARÁMETROS OPTIMIZADOS")
    print("=============================================")
    print(f"1. Módulo (q)        : {int(q):.4e}")
    print(f"   En potencias de 2 : 2^{gmpy2.log2(q):.2f}")
    print("---------------------------------------------")
    print(f"2. Ruido std (r/rho) : {r:.6f}")
    print(f"   (Debe ser >= 3.2) : {'CUMPLE' if r >= 3.2 else 'FALLA'}")
    print("=============================================")
    
    Aq = generate_Aq(q, coeffs)
    
    return Aq, tau, delta, m, N, coeffs, C_m, q, r

################################################################
#                           KeyGen                             #
################################################################

def KeyGen(Aq, N, r, q, rng, p):
    # Key Generation
    a = Aq.random_element()
    s_coeff, e_coeff= sample_discrete_gaussian_ZN(N, r, q, rng, num_samples=2)
    sk = Aq(s_coeff.tolist())
    e = Aq(e_coeff.tolist())
    b = (a*sk) + (p*e)

    pk = (a, b)
    # sk

    # Key hat Generation
    a_hat = Aq.random_element()
    b_hat = Aq.random_element()

    pk_hat = (a_hat, b_hat)
    
    return pk, sk, pk_hat

################################################################
#                    D_{\rho}^d  distribution                  #
################################################################

def sample_discrete_gaussian_ZN(N, s, q, rng, num_samples=1):
    """
    Muestrea desde D_{Z^N, s} con una única semilla global fija (SEED):
      1) x ~ R^N con densidad proporcional a exp(-pi * ||x||^2 / s^2)
      2) redondear componente a componente al entero más cercano
      3) reducir módulo q (devuelve representantes en 0..q-1)

    Parámetros:
      N           : dimensión del lattice Z^N
      s           : parámetro gaussiano
      q           : módulo para reducir
      num_samples : número de vectores a generar

    Retorno:
      array (num_samples, N) con valores en 0..q-1
    """
    sigma = s / np.sqrt(2.0 * np.pi)
    samples_cont = rng.normal(loc=0.0, scale=sigma, size=(num_samples, N))
    samples_rounded = np.rint(samples_cont).astype(np.int64)
    return np.mod(samples_rounded, q)

################################################################
#                       F_{p^k}  distribution                  #
################################################################

def sample_Fpk(p, k, n):
    # 1. Creamos el campo finito
    F = GF(p**k, 'a') 
    
    # 2. Generamos n elementos aleatorios
    # Usamos .vector() para convertir el polinomio (a^2+1) en lista [1, 0, 1...]
    return [F.random_element() for _ in range(n)]
