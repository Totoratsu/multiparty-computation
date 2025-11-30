import numpy as np
from functools import partial
from sage.all import Zmod, PolynomialRing, ZZ
from multiprocessing import Pool, cpu_count


################################################################
#                        Aq Generation                         #
################################################################

def parameters_worker(m, p, k, s):
    # IMPORTANTE: las importaciones deben hacerse dentro del worker para evitar
    # problemas de pickling al usar multiprocessing con Sage.
    from sage.all import Zmod, PolynomialRing, cyclotomic_polynomial, gcd, euler_phi

    # Queremos m tal que gcd(m, p) = 1 para que la teoría ciclotómica estándar funcione.
    # Cuando (p, m) = 1, la reducción de Φ_m modulo p factoriza en grados que dividen
    # el orden multiplicativo de p mod m — esto es clave para controlar los grados k'.
    if gcd(m, p) != 1:
        return None

    # Cota de seguridad que se usa en el paper
    phi_m = int(euler_phi(m))
    if phi_m < 14300:
        return None

    # Construimos el polinomio ciclotómico Φ_m en ℤ[x].
    Phi_ZZ = cyclotomic_polynomial(m)

    # Pasamos a R = 𝔽_p[x].
    R = PolynomialRing(Zmod(p), 'x')

    # Reducimos Φ_m mod p (coerción a R).
    Phi = R(Phi_ZZ)

    # Factorizamos Φ_m en 𝔽_p[x]; matemáticamente buscamos su descomposición en factores
    # irreducibles para obtener sus grados k'_i.
    fac = Phi.factor()

    # Construimos la lista de grados (teniendo en cuenta exponentes).
    # Esto produce los grados k'_i de los cuerpos finitos 𝔽_{p^{k'_i}}
    # que aparecen en la descomposición A ≅ ∏ 𝔽_{p^{k'_i}}.
    degrees = [f.degree() for (f, e) in fac for _ in range(e)]

    # Buscamos factores cuyo grado sea múltiplo de k, porque sólo así
    # 𝔽_{p^{k'_i}} contendrá una subextensión isomorfa a 𝔽_{p^k}.
    good_degs = [d for d in degrees if d % k == 0]

    # Si hay al menos s tales factores, entonces Φ_m satisface la condición para alojar
    # M = (𝔽_{p^k})^s dentro de A.
    if len(good_degs) >= s:
        # Para enviar datos entre procesos guardamos sólo objetos serializables.
        coeffs = [int(c) for c in Phi.list()]  # coeficientes de Φ_m mod p
        fac_serial = [(str(f), int(e)) for (f, e) in fac]
        return (m, coeffs, fac_serial, degrees, phi_m)

    return None

def calculate_m(p, k, s, min_m, max_m):
    # Preconfiguramos el worker para que reciba sólo m
    worker_partial = partial(parameters_worker, p=p, k=k, s=s)

    found = None
    pool = Pool(processes=cpu_count())
    try:
        # Recorremos los valores de m en paralelo.
        # Mat. hablando: probamos distintos índices ciclotómicos
        # buscando uno cuya reducción mod p produzca los grados deseados.
        for res in pool.imap_unordered(worker_partial, range(min_m, max_m + 1), chunksize=16):
            if res is not None:
                found = res
                break
    finally:
        pool.terminate()
        pool.join()

    if not found:
        print("No encontrado en rango; aumenta max_m o relaja condiciones.")
        return None

    m, coeffs, fac_serial, degrees, phi_m = found

    print("m =", m, "phi(m) =", phi_m, "degrees:", degrees)
    print("coeff:", coeffs)
    
    return m, phi_m, coeffs

def generate_Aq(q, coeffs):
# 1. Definir el anillo base Z_q (enteros modulo q)
    # q suele ser una potencia de 2 (ej. 2^64 o 2^128)
    R_q = Zmod(q)
    
    # 2. Definir el anillo de polinomios sobre Z_q
    PolyRing_q = PolynomialRing(R_q, 'x')
    
    # 3. Reconstruir el polinomio ciclotómico Phi_m usando los coeficientes
    # Los coeficientes vienen de Phi_m calculado en Z (o modulo p), 
    # pero aquí los interpretamos como elementos de Z_q.
    Phi = PolyRing_q(coeffs)

    # 4. Construir el anillo cociente A_q = Z_q[x] / (Phi_m)
    # Este es el espacio donde viven los textos cifrados.
    Aq = PolyRing_q.quotient(Phi, 'a')
    
    print(f"Construido A_q = {Aq}")
    
    return Aq

################################################################
#                         Others                               #
################################################################

def unique_prime_factors(n):
    """
    Calcula el conjunto de factores primos únicos de un número entero n.
    Ejemplo: 8 -> {2}, 12 -> {2, 3}, 7 -> {7}
    """
    if n <= 1:
        return set()

    factores = set()
    d = 2

    # Paso 1: Manejar el factor 2 (el único factor primo par)
    if n % d == 0:
        factores.add(d)
        while n % d == 0:
            n //= d

    # Paso 2: Manejar factores primos impares
    # Solo necesitamos verificar hasta la raíz cuadrada del número restante (n)
    d = 3
    limite = int(np.sqrt(n))
    while d <= limite:
        if n % d == 0:
            factores.add(d)
            while n % d == 0:
                n //= d
            # Recalcular el límite de la raíz cuadrada para el nuevo n reducido
            limite = int(np.sqrt(n))
        d += 2 # Saltar al siguiente impar (3, 5, 7, 9...)

    # Paso 3: Manejar el caso donde n es un primo grande
    # Si al final n es mayor que 1, el n restante es el último factor primo.
    if n > 1:
        factores.add(n)

    return factores

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

def get_Drhod(d):
    Drhod = []
    for _ in range(d):
        Drhod.append()