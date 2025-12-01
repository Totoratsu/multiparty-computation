# Clase Sage/Python para el anillo G := A^3_q con suma componente a componente
# y producto definido sólo cuando la 3ª componente de los dos factores es 0:
# (a0,a1,0)*(b0,b1,0) := (a0*b0, a1*b0 + a0*b1, - a1*b1)
#
# Uso previsto: pasar una instancia de A (por ejemplo A = R.quotient(Phi,'a') en Sage)
# y construir G = GRing(A). Luego crear elementos con G.element(...) o G(a0,a1,a2).
#
# Notas matemáticas en los comentarios:
# - A es un anillo (en nuestro contexto, una álgebra cociente sobre 𝔽_p).
# - G es el conjunto A^3 con estructura abeliana por suma componente a componente.
# - El producto está sólo definido cuando la tercera componente de ambos factores es 0;
#   la implementación lanzará excepción si se intenta multiplicar en otro caso.
# - Para facilitar la interoperabilidad aceptamos valores que sean convertibles a A
#   (por ejemplo enteros, polinomios en el anillo A, o elementos ya en A).

class GRing:
    def __init__(self, A):
        """
        A: el anillo/base (por ejemplo A = R.quotient(Phi,'a') en Sage).
        Guardamos A como parent ring para coerciones.
        """
        self.A = A
        # 3-tupla que representa el elemento neutro aditivo (0,0,0)
        self.zero = self.element(0, 0, 0)
        # definimos 1 multiplicativo natural como (1,0,0) si se quiere
        self.one = self.element(1, 0, 0)

    # fábrica / constructor conveniente
    def element(self, a0, a1, a2):
        """Crea un elemento de G a partir de tres componentes (coercibles a A)."""
        return GElement(a0, a1, a2, self)

    def __repr__(self):
        return "GRing(A=%s)" % (repr(self.A),)


class GElement:
    def __init__(self, a0, a1, a2, parent: GRing):
        """
        a0,a1,a2: componentes (serán coerced a parent.A)
        parent: instancia de GRing que contiene el anillo A
        """
        self.parent = parent
        A = parent.A
        # Coerción a A: esto permite pasar enteros, polinomios, o elementos de A.
        self.a0 = A(a0)
        self.a1 = A(a1)
        self.a2 = A(a2)

    # Representación legible
    def __repr__(self):
        return "GEl(%r, %r, %r)" % (self.a0, self.a1, self.a2)

    # Igualdad componente a componente (matemáticamente: igualdad en A^3)
    def __eq__(self, other):
        if not isinstance(other, GElement):
            return False
        if self.parent is not other.parent:
            # elementos de distintos anillos base no son comparables aquí
            return False
        return (self.a0 == other.a0) and (self.a1 == other.a1) and (self.a2 == other.a2)

    # Suma componente a componente (estructuralmente A^3 como grupo aditivo)
    def __add__(self, other):
        if not isinstance(other, GElement) or self.parent is not other.parent:
            raise TypeError("Suma sólo definida entre elementos del mismo GRing")
        return GElement(self.a0 + other.a0,
                        self.a1 + other.a1,
                        self.a2 + other.a2,
                        self.parent)

    def __radd__(self, other):
        # soporta sum builtin: 0 + x -> x
        if other == 0:
            return self
        return self.__add__(other)

    def __neg__(self):
        # inverso aditivo
        return GElement(-self.a0, -self.a1, -self.a2, self.parent)

    def __sub__(self, other):
        return self + (-other)

    # Producto según la fórmula dada:
    # sólo permitido si la tercera componente de ambos factores es 0.
    def __mul__(self, other):
        if not isinstance(other, GElement) or self.parent is not other.parent:
            raise TypeError("Producto sólo definido entre elementos del mismo GRing")
        A = self.parent.A
        zeroA = A(0)
        # la especificación: "el producto solo está definido cuando la tercera componente
        # de los dos vectores es 0"
        if (self.a2 != zeroA) or (other.a2 != zeroA):
            raise ValueError("Producto sólo definido si la tercera componente de ambos factores es 0")
        # calcular según la regla:
        # (a0,a1,0)*(b0,b1,0) = (a0*b0, a1*b0 + a0*b1, - a1*b1)
        c0 = self.a0 * other.a0
        c1 = self.a1 * other.a0 + self.a0 * other.a1
        c2 = -(self.a1 * other.a1)
        return GElement(c0, c1, c2, self.parent)

    # Support for left/right multiplication by scalars from A is NOT defined here
    # (no implementamos coerción automática de A->G). Si se desea, se puede añadir:
    def scalar_mul_left(self, scal):
        """Multiplicación por un escalar de A por la izquierda (escalar * elemento) componente a componente."""
        A = self.parent.A
        s = A(scal)
        return GElement(s * self.a0, s * self.a1, s * self.a2, self.parent)

    def scalar_mul_right(self, scal):
        """Multiplicación por un escalar de A por la derecha (elemento * escalar) componente a componente."""
        A = self.parent.A
        s = A(scal)
        return GElement(self.a0 * s, self.a1 * s, self.a2 * s, self.parent)

    # Opcional: norma simple/validación
    def is_mul_defined_with(self, other):
        """Devuelve True si el producto con `other` está definido (terceras componentes = 0)."""
        A = self.parent.A
        return (self.a2 == A(0)) and (other.a2 == A(0))


# ---------------------------
# Ejemplo de uso (en Sage)
# ---------------------------
# Asumiendo que ya tienes A creado en el entorno Sage, por ejemplo:
#   R = PolynomialRing(Zmod(p), 'x')
#   Phi = R(cyclotomic_polynomial(m))
#   A = R.quotient(Phi, 'a')
#
# entonces:
#
#   G = GRing(A)
#   x = G.element(1, 2, 0)        # crea (1,2,0) con coerción a A
#   y = G.element(3, 5, 0)        # crea (3,5,0)
#   z = x * y                     # aplica la regla de producto; funciona porque 3ª comp = 0
#   s = x + y                     # suma componente a componente
#   print(z)                      # muestra resultado en términos de elementos de A
#
# Si intentas multiplicar elementos cuya 3ª componente no es 0:
#   bad = G.element(1,2,7)
#   x * bad   # -> ValueError: Producto sólo definido si la tercera componente de ambos factores es 0
#
# Eso es todo: la clase encapsula la estructura aditiva de A^3 y la regla de producto que pediste.
