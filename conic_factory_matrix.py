import numpy as np
import matplotlib.pyplot as plt
import sympy
from fractions import Fraction
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
            return Circle(equation, conic.coefficients)
        elif conic_type == "Ellipse":
            return Ellipse(equation, conic.coefficients)
        elif conic_type == "Parabola":
            return Parabola(equation, conic.coefficients)
        elif conic_type == "Hyperbola":
            return Hyperbola(equation, conic.coefficients)
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
        B = self.coeff.get((1, 1), 0) * 2  # This is because the coefficient of xy in the equation is B/2
        C = self.coeff.get((0, 2), 0)
        D = self.coeff.get((1, 0), 0) * 2  # This is because the coefficient of x in the equation is D/2
        E = self.coeff.get((0, 1), 0) * 2  # This is because the coefficient of y in the equation is E/2
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

    def _to_general_form_conic(self, equation: str):
        """
        Takes the string, works out the LCM, collects like terms
        then reduces the resultant as far as possible.
        :return: General form equation of a conic as string and sympy expression
        """
        x, y = sympy.symbols('x y')
        equation_str = sympy.sympify(equation)
        equation_str = sympy.expand(equation_str)
        equation_str = self._remove_fraction(equation_str)
        equation_str = sympy.collect(equation_str, (x, y))
        formula = sympy.simplify(equation_str)
        formula_str = str(formula)

        return formula_str, formula

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
    def __init__(self, equation, coefficients):
        super().__init__(equation)

    def draw(self):
        pass  # Implement plotting for ellipse


class Ellipse(Conic):
    def __init__(self, equation, coefficients):
        super().__init__(equation)
        # Compute any properties unique to ellipses here

    def draw(self):
        pass  # Implement plotting for ellipse


class Circle(Conic):
    def __init__(self, equation, coefficients):
        super().__init__(equation)
        # Compute any properties unique to circles here

    def draw(self):
        pass  # Implement plotting for circle


class Hyperbola(Conic):
    def __init__(self, equation, coefficients):
        super().__init__(equation)
        # Compute any properties unique to hyperbolas here

    def draw(self):
        pass  # Implement plotting for hyperbola



