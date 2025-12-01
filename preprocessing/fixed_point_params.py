import numpy as np
import gmpy2
from gmpy2 import mpfr

# -------------------------------------------------------
# Cambios mínimos: usar gmpy2 solo donde hay potencias/log/ceil grandes
# -------------------------------------------------------

# =========================================================
# 2. FUNCIONES DE RESTRICCIÓN (SOLO gmpy2 donde hace falta)
# =========================================================

def constraint_q(r_i, N, C_m, p, tau, c_sec, n, sec):
    """
    Calcula el límite inferior necesario para q dado un ruido r_i.
    Cambios mínimos: se reemplaza solo el cálculo de 2^sec por una operación
    exacta (shift) cuando sec es entero, evitando np.power(2.0, sec).
    """
    # Término cuadrático del ruido (r^2 * N^2)
    noise_quad = np.power(r_i * N, 2)

    # Fórmula de Y (Ruido acumulado intermedio)
    term_1 = 4 * C_m * noise_quad
    term_2 = 2 * np.sqrt(N) * r_i
    term_3 = 4 * C_m * noise_quad

    Y = tau + p * (term_1 + term_2 + term_3)

    # Fórmula de Z (Límite de ruido final B_final)
    factor_comun = n * N * c_sec * Y
    Z = C_m * np.power(factor_comun, 2) + n * c_sec * Y

    # Límite final de q: q > 2 * Z * (1 + 2^sec)
    # --- Cambio mínimo: calcular 2^sec de forma exacta usando gmpy2 / shift
    if isinstance(sec, int):
        pow2sec = 1 << sec   # uso de shift: exacto y rápido (Python int arbitrary precision)
    else:
        # caso no habitual si sec no es entero: usar gmpy2 para potencia en punto flotante
        pow2sec = float(gmpy2.mpfr(mpfr(2))**mpfr(sec))

    constant_q = 2 * (1 + pow2sec)
    q_limit = constant_q * Z

    return q_limit


def constraint_r(q_i, N, delta):
    """
    Calcula el límite inferior para r_std (rho) dado un módulo q_i.
    Basado en la dificultad del problema LWE.
    """
    # Prevenir errores si q_i es muy pequeño
    if q_i <= 1:
        return 3.2

    # Convertir a mpfr solo para las operaciones log/expcriticas
    q_m = mpfr(q_i)
    delta_m = mpfr(delta)
    N_m = mpfr(N)  # N es entero, lo pasamos a mpfr solo para multiplicar con mpfr

    # log2 usando gmpy2 (más estable para grandes)
    log_q = gmpy2.log2(q_m)
    log_delta = gmpy2.log2(delta_m)

    # t' = sqrt( N * log2(q) / log2(delta) )
    t_prime = gmpy2.sqrt(N_m * log_q / log_delta)

    # exponente = 1 - N/t'
    exponente = mpfr(1) - (N_m / t_prime)

    # var = 1.5 * delta^(-t') * q^(1 - N/t')
    term_delta = gmpy2.mpfr(delta_m)**(-t_prime)   # delta ** (-t_prime)
    term_q = gmpy2.mpfr(q_m)**exponente          # q_i ** exponente

    var_mp = mpfr('1.5') * term_delta * term_q

    # Devolver como float (comportamiento igual al original)
    var = float(var_mp)
    r_limit = max(3.2, var)

    return r_limit


def fixed_point_qr(r_initial, N, C_m, p, tau, c_sec, n, sec, delta, max_iter=200, tolerance=1e-4):
    """
    Encuentra el punto de equilibrio (Fixed Point) entre las restricciones de q y r.
    Cambios mínimos: cálculo de la siguiente potencia de 2 usando gmpy2.log2/ceil o shift.
    """
    r_i = r_initial
    q_i = 0.0

    print(f"--- Iniciando Iteración (N={N}, sec={sec}) ---")

    for i in range(max_iter):
        # A. Calcular cota inferior de q necesaria para la CORRECCIÓN (dado r actual)
        q_min_needed = constraint_q(r_i, N, C_m, p, tau, c_sec, n, sec)

        # Ajustamos q al siguiente valor "válido" (generalmente potencia de 2 o primo cercano)
        # Aquí usamos la siguiente potencia de 2:
        # --- Cambio mínimo: usar gmpy2 para log2 y ceil cuando q_min_needed grande ---
        e_mp = gmpy2.ceil(gmpy2.log2(mpfr(q_min_needed)))
        e = int(e_mp)
        q_i = gmpy2.mpz(2)**e if e < 1024 else (1 << e)  # si e grande, usar shift para exactitud

        # B. Calcular cota inferior de r necesaria para la SEGURIDAD (dado q actual)
        r_new = constraint_r(q_i, N, delta)

        # C. Chequeo de Convergencia
        diff = np.abs(r_new - r_i)
        if diff < tolerance:
            print(f"✅ Convergencia en iteración {i+1}")
            return q_i, r_new

        # Actualizar r para la siguiente iteración
        r_i = r_new

        if i % 10 == 0:
            print(f"Iter {i}: r_std={r_i:.4f} => q_necesario=2^{np.log2(q_i):.2f}")

    print("⚠️ No se alcanzó convergencia estable (revisar parámetros N o delta).")
    return q_i, r_i
