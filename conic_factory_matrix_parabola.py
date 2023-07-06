import numpy as np
import matplotlib.pyplot as plt
import sympy
from fractions import Fraction
from math import sin, cos
from IPython.display import display
from poly_dictionary import decompose_polynomial


class Conic:
    def __init__(self, equation: str):
        self.original_input = equation
        self.equation_str, self.equation = self._to_general_form_conic(equation)
        self.coefficients = decompose_polynomial(self.equation_str)
        self.coefficients = self._remove_z_coordinate(self.coefficients)
        self.coeff_matrix = self.create_coeff_matrix()
        self.type = self.classify()
        self.expression = None

    @classmethod
    def create(cls, equation: str):
        # Create an instance of Conic
        conic = cls(equation)

        # Based on the classify method, create appropriate subclass
        conic_type = conic.classify()
        if conic_type == "Circle":
            return Circle(equation)
        elif conic_type == "Ellipse":
            return Ellipse(equation)
        elif conic_type == "Parabola":
            return Parabola(equation)
        elif conic_type == "Hyperbola":
            return Hyperbola(equation)
        else:
            raise ValueError('Could not determine the type of the conic section')

    @staticmethod
    def _remove_z_coordinate(coefficients):
        """
        Removes any terms with a z coordinate from the coefficient dictionary.
        :param : Original dictionary of form {(x, y, z): coefficient, ... }
        :return: New dictionary of form {(x, y): coefficient, ... }
        """
        return {(x, y): coeff for (x, y, z), coeff in coefficients.items() if z == 0}

    @staticmethod
    def _remove_fraction(equation_str):
        """
        Removes any fractional terms from the equation.
        """
        denominators = [sympy.fraction(term)[1] for term in equation_str.as_ordered_terms()]

        if all(denom == 1 for denom in denominators):
            return str(equation_str)

        lcm_value = sympy.lcm(denominators)
        equation_no_frac = equation_str * lcm_value

        return equation_no_frac

    @property
    def coeff(self):
        return self.coefficients

    def __repr__(self):
        return f"Original : {self.original_input}"

    def __str__(self):
        """
        Return the equation in user-friendly syntax
        """
        equation = sympy.printing.sstr(self.equation)
        return equation.replace('**', '^').replace('*', '')

    def create_coeff_matrix(self):
        """
        Create the matrix of coefficients, as you suggested.
        We can then use this matrix in all subsequent calculations.
        """
        A = self.coeff.get((2, 0), 0)
        B = self.coeff.get((1, 1), 0)   # This is because the coefficient of xy in the equation is B/2
        C = self.coeff.get((0, 2), 0)
        D = self.coeff.get((1, 0), 0)   # This is because the coefficient of x in the equation is D/2
        E = self.coeff.get((0, 1), 0)   # This is because the coefficient of y in the equation is E/2
        F = self.coeff.get((0, 0), 0)

        matrix = np.array([[A, B / 2, D / 2], [B / 2, C, E / 2], [D / 2, E / 2, F]])
        return matrix

    def print_matrix(self):
        print(self.coeff_matrix)

    def classify(self):
        """
        Classify the conic section based on its coefficients.
        Using the determinant of the 2x2 sub-matrix (delta) and the determinant of the 3x3 matrix (Delta).
        """
        # Calculate the determinant of the 2x2 sub-matrix
        delta = np.linalg.det(self.coeff_matrix[:2, :2])

        # Calculate trace of the 2x2 sub-matrix
        tau = np.trace(self.coeff_matrix[:2, :2])

        # Calculate the determinant of the 3x3 matrix. Defined but unused... could be useful
        Delta = np.linalg.det(self.coeff_matrix)

        A = self.coeff_matrix[0, 0]
        B = self.coeff_matrix[0, 1] * 2
        C = self.coeff_matrix[1, 1]

        if delta == 0:
            return "Parabola"
        elif delta < 0:
            return "Hyperbola"
        elif delta > 0:
            if A == C and B == 0:
                return "Circle"
            else:
                return "Ellipse"
        else:
            return "Unknown"

    def _to_general_form_conic(self, equation: str):
        """
        Takes the string, works out the LCM, collects like terms
        then reduces the resultant as far as possible.
        :return: General form equation of a conic as string and sympy expression
        """
        x, y = sympy.symbols('x y')
        equation_str = sympy.sympify(equation)

        # we have in a form y = ax^2 + dx + f or x = cy^2 + ey + f
        # WHAT HAPPENS, when we set all to zero?
        # subtract either x or y from either side!!

        # Check if "y" term exists
        if 'y' not in str(equation_str):
            equation_str -= y  # subtracting y term - quadratics of the form x^2 + 5x + 9

        # Check if "x" term exists
        if 'x' not in str(equation_str):
            equation_str -= x  # subtracting x term - quadratics of form y^2 + 5y + 9

        equation_str = sympy.expand(equation_str)
        equation_str = self._remove_fraction(equation_str)
        equation_str = sympy.collect(equation_str, (x, y))

        formula = sympy.simplify(equation_str)
        formula_str = str(formula)

        return formula_str, formula

    def _convert_to_sympy_expression(self, rational=False, force=True):
        """
        Converts the current polynomial expression into a sympy expression.
        """
        if self.expression is None or force:
            # Create the symbols x, y
            x, y = sympy.symbols('x y')

        if rational:
            polynomial = sum(
                sympy.Rational(str(coeff)) * x ** powers[0] * y ** powers[1] for powers, coeff in self.coeff.items())
        else:
            polynomial = sum(
                coeff * x ** powers[0] * y ** powers[1] for powers, coeff in self.coeff.items())

        # Create a sympy.Poly object
        self.expression = sympy.Poly(polynomial, (x, y))

        return self.expression

    def display(self):
        """
        Displays the sympy expression of the polynomial using IPython's display system.
        """
        self._convert_to_sympy_expression()
        display(self.expression)

    def save_as_sympy(self, rational=False, return_expr=False):
        """
        Saves the polynomial as a sympy expression and returns it.
        """
        self._convert_to_sympy_expression(rational, force=True)
        if return_expr:
            return self.expression.as_expr()  # Return as a sympy expression
        else:
            return self.expression  # Return as a sympy.Poly object

    def plot(self):
        A, B, C, D, E, F = self.coeff_matrix[0, 0], self.coeff_matrix[0, 1] * 2, self.coeff_matrix[1, 1], \
                           self.coeff_matrix[0, 2] * 2, self.coeff_matrix[1, 2] * 2, self.coeff_matrix[2, 2]

        x = np.linspace(-10, 10, 400)
        y = np.linspace(-10, 10, 400)
        x, y = np.meshgrid(x, y)

        plt.contour(x, y, (A * x ** 2 + B * x * y + C * y ** 2 + D * x + E * y + F), [1], colors='r')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.axhline(0, color='gray', linewidth=0.5)
        plt.axvline(0, color='gray', linewidth=0.5)
        plt.title(f"${self.__str__()}$")
        plt.show()

    def compute_transformation(self):
        """
        This will probably be subclass method
        :return: the affine transformation mapping to standard form including the matrix required
        """
        # Compute the transformation matrix
        pass

    def to_standard_form(self):
        """
        Also subclass method
        :return: equation of conic in standard form
        """
        # Apply the transformation matrix to the equation and convert it to standard form
        pass


class Parabola(Conic):
    def __init__(self, equation):
        super().__init__(equation)
        self.A = self.coeff_matrix[0, 0]      # x^2 term
        self.B = self.coeff_matrix[0, 1] * 2  # Coefficient of xy
        self.C = self.coeff_matrix[1, 1]      # y^2 term
        self.D = self.coeff_matrix[0, 2] * 2  # Coefficient of x
        self.E = self.coeff_matrix[1, 2] * 2  # Coefficient of y
        self.F = self.coeff_matrix[2, 2]      # Constant
        self.orientation = self.get_orientation()
        self.vertex = self.compute_vertex()
        self.complete_sq_coeff = self.complete_sq_coeff()

    def get_orientation(self):
        # Determine the orientation based on the values of A, B, and C
        if self.B == 0:
            if self.A > 0 and self.C == 0:
                return "vertical", "positive"
            elif self.A == 0 and self.C > 0:
                return "horizontal", "negative"
            elif self.A < 0 and self.C == 0:
                return "vertical", "negative"
            elif self.A == 0 and self.C < 0:
                return "horizontal", "positive"
        else:
            # For a rotated parabola, return 'R' for now.
            # You will need to compute the angle of rotation here.
            theta = ...  # compute theta based on the values of A, B, and C
            return "rotated", theta

    def rotate(self):
        if self.orientation[0] == "vertical":
            if self.orientation[1] == "positive":
                theta = 3 * np.pi / 4
            else:  # "negative"
                theta = np.pi / 2
        elif self.orientation[0] == "horizontal":
            if self.orientation[1] == "positive":
                theta = 0
            else:  # "negative"
                theta = np.pi
        else:  # "rotated"
            theta = self.orientation[1]  # theta was stored in orientation tuple

        rotation_matrix = np.array([[np.cos(theta), np.sin(theta)],
                                    [-np.sin(theta), np.cos(theta)]])
        rotation = self.matrix().T @ rotation_matrix @ self.matrix
        return rotation

    # TODO: Can we just rotate first based on orientation? then we just need to map to the origin...
    # TODO: Need some way of storing the new matrix based on the rotation.
    # should directly compare this to the instantiation equation.

    def compute_vertex(self):
        if self.orientation[0] == 'vertical':  # Parabola opens up or down
            a = self.A / -self.E
            d = self.D / -self.E
            f = self.F / -self.E
            h = sympy.Rational(-d / (2 * a)).limit_denominator(100000)
            k = sympy.Rational(a * h ** 2 + d * h + f).limit_denominator(100000)
        else:  # Parabola opens to the right or left
            c = self.C / -self.D
            e = self.E / -self.D
            f = self.F / -self.D
            k = sympy.Rational(-e / (2 * c)).limit_denominator(100000)
            h = sympy.Rational(c * k ** 2 + e * k + f).limit_denominator(100000)
        return h, k

    def complete_sq_coeff(self):
        # Parabola opens upwards or downwards
        if self.B == 0:
            a = sympy.Rational(-self.E).limit_denominator(100000)
        else:  # Parabola opens to the right or left
            a = sympy.Rational(-self.D).limit_denominator(100000)
        return a

    def complete_square(self):
        """
        This method irrelevant really as we can compute the vertex.
        Why do we need this...? Not sure,
        nice for humans?
        :return:
        """
        h, k = self.vertex
        a = self.complete_sq_coeff

        h_sign = '-' if h < 0 else '+'
        k_sign = '-' if k < 0 else '+'
        h, k = abs(h), abs(k)

        if self.orientation[0] == 'vertical':  # Parabola opens upwards or downwards
            return f"y = (1 / {a}) * (x {h_sign} {h})^2 {k_sign} {k}".replace(" - - ", " + ").replace(" + - ", " - ")
        else:  # Parabola opens to the right or left
            return f"x = (1 / {a}) * (y {k_sign} {k})^2 {h_sign} {h}".replace(" - - ", " + ").replace(" + - ", " - ")

    def plot(self):
        h, k = self.vertex
        margin = 5  # adjust this to your needs

        if self.orientation[0] == 'vertical':  # Parabola opens up or down
            x = np.linspace(float(h) - margin, float(h) + margin, 400)
            y = (self.A * x ** 2 + self.D * x + self.F) / -self.E
        else:  # Parabola opens to the right or left
            y = np.linspace(float(k) - margin, float(k) + margin, 400)
            x = (self.C * y ** 2 + self.E * y + self.F) / -self.D

        plt.plot(x, y, color='r')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.axhline(0, color='gray', linewidth=0.5)
        plt.axvline(0, color='gray', linewidth=0.5)
        plt.plot(h, k, 'bo')  # plot the vertex as a blue dot
        plt.annotate(f'({h}, {k})', (h, k), textcoords="offset points", xytext=(-10, -10), ha='center')  # Annotate vertex
        plt.title(f"Completed Square:\n${self.complete_square()}$")
        plt.show()


class Ellipse(Conic):
    def __init__(self, equation):
        super().__init__(equation)
        # Compute any properties unique to ellipses here

    def draw(self):
        pass  # Implement plotting for ellipse


class Circle(Conic):
    def __init__(self, equation):
        super().__init__(equation)
        # Compute any properties unique to circles here

    def draw(self):
        pass  # Implement plotting for circle


class Hyperbola(Conic):
    def __init__(self, equation):
        super().__init__(equation)
        # Compute any properties unique to hyperbolas here

    def draw(self):
        pass  # Implement plotting for hyperbola



