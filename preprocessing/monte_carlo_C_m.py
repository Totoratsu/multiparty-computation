import numpy as np
import multiprocessing
from functools import partial
from sage.all import ZZ, PolynomialRing, cyclotomic_polynomial, euler_phi

# --- 1. FUNCIÓN WORKER (Cálculo con Norma L2) ---
def worker_montecarlo_l2(batch_size, m):
    """
    Calcula el factor de expansión C_m usando la norma L2 (Euclidiana), 
    que es la estándar en el análisis de ruido en LWE.
    """
    # Importaciones requeridas para el proceso hijo
    from sage.all import ZZ, PolynomialRing, cyclotomic_polynomial, euler_phi
    import numpy as np
    
    try:
        R = PolynomialRing(ZZ, 'x')
        Phi = cyclotomic_polynomial(m)
        n = euler_phi(m)
        
        local_max = 0.0
        
        for _ in range(batch_size):
            # Generamos polinomios con coeficientes {-1, 0, 1}
            a = R.random_element(degree=n-1, x=-1, y=2)
            b = R.random_element(degree=n-1, x=-1, y=2)
            
            if a.is_zero() or b.is_zero():
                continue
                
            # Operación: Multiplicación y Reducción
            res = (a * b) % Phi
            
            # 2. CALCULAR NORMAS L2 (Euclidianas)
            
            # Conversión a arrays de numpy para calcular normas más rápido
            coeffs_a = np.array([int(c) for c in a.list()])
            coeffs_b = np.array([int(c) for c in b.list()])
            coeffs_res = np.array([int(c) for c in res.list()])
            
            # Norma L2: sqrt(suma de cuadrados)
            norm_a_l2 = np.linalg.norm(coeffs_a, 2)
            norm_b_l2 = np.linalg.norm(coeffs_b, 2)
            norm_res_l2 = np.linalg.norm(coeffs_res, 2)
            
            # Prevenir división por cero
            if norm_a_l2 == 0 or norm_b_l2 == 0:
                continue
                
            # 3. Ratio de expansión L2
            # C_m = ||res||_2 / (||a||_2 * ||b||_2)
            ratio = norm_res_l2 / (norm_a_l2 * norm_b_l2)
            
            if ratio > local_max:
                local_max = ratio
                
        return local_max

    except Exception as e:
        print(f"Error en worker para m={m}: {e}")
        return 0.0

def calculate_Cm_parallel(m, total_trials=1000):
    # m es asumido como entero
    num_cores = multiprocessing.cpu_count()
    
    batch_size = total_trials // num_cores
    tasks = [batch_size] * num_cores
    remainder = total_trials % num_cores
    if remainder > 0:
        tasks[-1] += remainder
        
    print(f"--- Iniciando Monte Carlo Paralelo (L2 Norm) para m={m} ---")
    print(f"--- Repartiendo {total_trials} pruebas en {num_cores} núcleos ---")

    # USAMOS EL NUEVO WORKER
    worker_with_args = partial(worker_montecarlo_l2, m=m)
    
    with multiprocessing.Pool(processes=num_cores) as pool:
        results = pool.map(worker_with_args, tasks)
    
    global_max = max(results)
    return global_max

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

    return list(factores)