Fase Online MPC + Reshare — Documentación Técnica
Este documento resume toda la implementación realizada de la fase online del protocolo SPDZ, incluyendo:

* Representación de authenticated shares (AuthShare)
* Subrutinas principales:Input, Multiply, Output
* Protocolo Reshare (opciones preserve-secret y delta-refresh)
* Breve descripción del flujo completo
* Interfaz de red simulada
* Cómo conectar esto con el preprocesamiento existente

1. Arquitectura General
La fase online del protocolo MPC (tipo SPDZ) utiliza dos ingredientes:

Shares autenticados
Cada valor secreto x se almacena como un objeto:


<x> = (shares[i], macs[i])   para i = 0..n-1
con: sum(shares[i]) = x + δ
     sum(macs[i])   = α · (x + δ)



donde α es la clave global secreta del MAC y δ es un offset aleatorio.
2. Fase Online
Usa el preprocessing para ejecutar:

* Input de valores
* Multiplicación
* Output (apertura segura)
* Reshare (recableado seguro de shares)

 2. AuthShare: representación de secretos autenticados
 class AuthShare:
    def __init__(self, F, delta, shares, macs):
        self.F = F
        self.delta = F(delta)
        self.shares = [F(s) for s in shares]
        self.macs  = [F(m) for m in macs]

    def value(self):
        return sum(self.shares)

    def mac_sum(self):
        return sum(self.macs)

    def n_players(self):
        return len(self.shares)
Este objeto encapsula:
* las shares
* los MAC shares
* el delta
* la aritmética en el campo GF(p)

3. Input Phase — Protocolo
El objetivo:
Un participante quiere introducir un valor público x, pero que quede en forma autenticada <x>.
Flujo (SPDZ clásico):

Cliente y jugadores tienen un <r> del preprocessing

Cliente abre ε = x - r
Jugadores computan: <x> = <r> + ε

Implementación 
def input_phase(F, alpha, x_input, R):
    r = R.value()

    eps = F(x_input - r)

    new_shares = [R.shares[i] + eps for i in range(n)]
    new_macs   = [R.macs[i] for i in range(n)]
    new_macs[-1] += eps * alpha

    return AuthShare(F, R.delta, new_shares, new_macs), eps
Esto genera un <x> correcto y autenticado. 

4. Multiply Phase — Multiplicación con triple de Beaver
Usamos un triple <a>, <b>, <c> con c = a·b.

Protocolo:
e = x - a
d = y - b
Abrir e y d
<x*y> = <c> + e<b> + d<a> + e·d
Implementación:
def multiply_phase(F, X, Y, triple):
    A,B,C = triple

    e = X.value() - A.value()
    d = Y.value() - B.value()

    Z_sh = [C.shares[i] + e*B.shares[i] + d*A.shares[i] for i in range(n)]
    Z_mac = [C.macs[i] + e*B.macs[i] + d*A.macs[i] for i in range(n)]

    Z = AuthShare(F, C.delta, Z_sh, Z_mac)
    return Z, e, d
5. Output Phase — Apertura segura del resultado
Se abre <z> verificando:
sum(macs_i) == α * sum(shares_i)
Si la verificación MAC pasa, el valor es correcto.
def output_phase(F, alpha, opened_list, Z):
    z_val = Z.value()
    mac_sum = Z.mac_sum()

    if mac_sum != alpha * z_val:
        return False, None

    return True, z_val
 6. Protocolo RESHARE — Realeatorizar shares
El reshare permite:

Reexpresar <x> bajo nuevas shares

Mantener consistencia del MAC

Opcional:

preserve-secret: garantiza que el valor no cambie

delta-refresh: produce un delta nuevo más seguro

6.1 Implementación (versión práctica)
def reshare_practical(F, X_auth, R_auth, alpha, network, preserve_secret=True):
    n = X_auth.n_players()
    network.start_round()

    d_i = [X_auth.shares[i] - R_auth.shares[i] for i in range(n)]
    for i in range(n):
        network.broadcast(i, d_i[i])

    msgs = network.collect()
    eps = F(sum(v for (_,v) in msgs))

    r = R_auth.value()
    x = r + eps

    if preserve_secret:
        delta_new = F(0)
        target = x
    else:
        delta_new = F.random_element()
        target = x + delta_new

    new_shares = [F.random_element() for _ in range(n-1)]
    new_shares.append(target - sum(new_shares))

    mac_target = alpha * target
    new_macs = [F.random_element() for _ in range(n-1)]
    new_macs.append(mac_target - sum(new_macs))

    network.end_round()

    return AuthShare(F, delta_new, new_shares, new_macs), {"eps": int(eps)}
    | Caso                    | Uso                                                     |
| ----------------------- | ------------------------------------------------------- |
| `preserve_secret=True`  | Para refrescar shares sin alterar el valor              |
| `preserve_secret=False` | Para ocultar aún más el valor introduciendo nuevo delta |

7. Red simulada (SimNetwork)
Proporciona una abstracción que permite simular broadcast y recepción.
class SimNetwork:
    def __init__(self, n): self.n = n
    def start_round(self): self.buffer=[]
    def broadcast(self, sender, value): self.buffer.append((sender,value))
    def collect(self): return self.buffer[:]
    def end_round(self): pass

8. Flujo completo integrado (online + reshare)
# 1. Reshare para preparar estado
X2,_ = reshare_practical(F, X, Rmask, alpha, net, preserve_secret=True)

# 2. Input de usuario
X_input, eps1 = input_phase(F, alpha, user_value, Rmask)

# 3. Multiplicación
Z, e, d = multiply_phase(F, X_input, Y_input, triple)

# 4. Output
ok, result = output_phase(F, alpha, [X_input, Y_input], Z)
 9. Conexión con el Preprocesamiento
Tu preprocesamiento produce:

<r> máscaras

<a>, <b>, <c> triples

alpha

parámetros del campo

anillo para cifrado homomórfico

La fase online solo requiere:
* el campo GF(p)
* las máscaras <r>
* los triples
* la clave alpha

Por eso ya podemos enlazar completamente ambas fases.

10. Verificación
Utilidad auxiliar:
def verify_authshare_consistency(A, alpha):
    return int(A.mac_sum()) == int(alpha * A.value())
