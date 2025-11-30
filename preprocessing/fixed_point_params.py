import numpy as np


# =========================================================
# 2. FUNCIONES DE RESTRICCIÓN (Corregidas)
# =========================================================

def constraint_q(r_i, N, C_m, p, tau, c_sec, n, sec):
    """
    Calcula el límite inferior necesario para q dado un ruido r_i.
    Basado en el crecimiento del ruido en el circuito aritmético.
    """
    # CORRECCIÓN: En tu código original tenías '4*C_m**np.pow(...)'.
    # Asumí que era una multiplicación (4 * C_m * ...) y no una exponenciación enorme.
    # También unifiqué los términos semejantes si la fórmula lo permite.
    
    # Término cuadrático del ruido (r^2 * N^2)
    noise_quad = np.power(r_i * N, 2)
    
    # Fórmula de Y (Ruido acumulado intermedio)
    # Interpretación: tau + p * (términos de ruido de multiplicación)
    term_1 = 4 * C_m * noise_quad
    term_2 = 2 * np.sqrt(N) * r_i
    term_3 = 4 * C_m * noise_quad # Este término aparecía repetido en tu fórmula original
    
    Y = tau + p * (term_1 + term_2 + term_3)
    
    # Fórmula de Z (Límite de ruido final B_final)
    # Z = C_m * (n * N * c_sec * Y)^2 + n * c_sec * Y
    factor_comun = n * N * c_sec * Y
    Z = C_m * np.power(factor_comun, 2) + (n * c_sec * Y)
    
    # Límite final de q: q > 2 * Z * (1 + 2^sec)
    constant_q = 2 * (1 + np.power(2.0, sec))
    q_limit = constant_q * Z
    
    return q_limit

def constraint_r(q_i, N, delta):
    """
    Calcula el límite inferior para r_std (rho) dado un módulo q_i.
    Basado en la dificultad del problema LWE.
    """
    # Prevenir errores matemáticos si q_i es muy pequeño (no debería ocurrir)
    if q_i <= 1: return 3.2
    
    log_q = np.log2(q_i)
    log_delta = np.log2(delta)
    
    # 1. Cálculo del exponente t' (t_prime)
    # t' = sqrt( N * log2(q) / log2(delta) )
    t_prime = np.sqrt(N * log_q / log_delta)
    
    # Verificación de seguridad: N debe ser <= t' para que la seguridad tenga sentido en esta fórmula?
    # Revisando tu fórmula: exponente = 1 - N/t'. Si N > t', el exponente es negativo.
    # Para q grande, exponente negativo hace que 'var' sea muy pequeño, dominando 3.2.
    # Esto es el comportamiento esperado: si q es gigante, r puede ser pequeño.
    
    # 2. Término de Seguridad
    exponente = 1.0 - (N / t_prime)
    
    # var = 1.5 * delta^(-t') * q^(1 - N/t')
    term_delta = np.power(delta, -t_prime)
    term_q = np.power(q_i, exponente)
    
    var = 1.5 * term_delta * term_q
    
    # 3. Límite inferior (Max entre 3.2 y el requisito LWE)
    r_limit = np.max([3.2, var])
    
    return r_limit

def fixed_point_qr(r_initial, N, C_m, p, tau, c_sec, n, sec, delta, max_iter=200, tolerance=1e-4):
    """
    Encuentra el punto de equilibrio (Fixed Point) entre las restricciones de q y r.
    """
    r_i = r_initial
    q_i = 0.0
    
    print(f"--- Iniciando Iteración (N={N}, sec={sec}) ---")
    
    for i in range(max_iter):
        # A. Calcular cota inferior de q necesaria para la CORRECCIÓN (dado r actual)
        q_min_needed = constraint_q(r_i, N, C_m, p, tau, c_sec, n, sec)
        
        # Ajustamos q al siguiente valor "válido" (generalmente potencia de 2 o primo cercano)
        # Aquí usamos la siguiente potencia de 2 para simplificar.
        q_i = np.power(2.0, np.ceil(np.log2(q_min_needed)))
        
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
