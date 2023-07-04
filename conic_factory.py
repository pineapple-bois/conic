import numpy as np
import matplotlib.pyplot as plt
import sympy
from fractions import Fraction
from poly_dictionary import decompose_polynomial


class Conic:
    def __init__(self, equation):
        self.original_input = equation
        self.equation = self._to_general_form_conic(equation)
        self.coefficients = decompose_polynomial(self.equation)
        self.coefficients = self._remove_z_coordinate(self.coefficients)
        self.type = self.classify()  # Determine the type of conic section based on coefficients

    @classmethod
    def create(cls, equation):
        # Create an instance of Conic
        conic = cls(equation)

        # Unpack coefficients for the general form
        A = conic.coefficients.get((2, 0), 0)
        B = conic.coefficients.get((1, 1), 0)
        C = conic.coefficients.get((0, 2), 0)
        D = conic.coefficients.get((1, 0), 0)
        E = conic.coefficients.get((0, 1), 0)
        F = conic.coefficients.get((0, 0), 0)

        discriminant = B**2 - 4*A*C

        # Based on the discriminant, create appropriate subclass
        if discriminant < 0 or (A == C and B == 0):
            # If A = C and B = 0, it's a circle regardless of the discriminant
            return Circle(equation, conic.coefficients) if A == C and B == 0 \
                            else Ellipse(equation, conic.coefficients)
        elif discriminant == 0:
            return Parabola(equation, conic.coefficients)
        elif discriminant > 0:
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

    @property
    def coeff(self):
        return self.coefficients

    def _to_general_form_conic(self, equation: str):
        """
        takes the string, works out the LCM, collects like terms
        then reduces the resultant as far as possible.
        :return: General form equation of a conic as string
        """
        # Create symbols
        x, y = sympy.symbols('x y')
        try:
            # Convert the input string to a sympy expression
            equation_str = sympy.sympify(equation)
        except sympy.SympifyError:
            raise ValueError("Invalid input. Please provide a valid sympy expression.")

        # Expand the expression
        equation_str = sympy.expand(equation_str)

        # Collect the denominators from the equation
        denominators = [sympy.fraction(term)[1] for term in equation_str.as_ordered_terms()]

        # If the equation already has no fractions, return as is
        if all(denom == 1 for denom in denominators):
            return str(equation_str)  # Changed here

        # Compute the least common multiple of the denominators
        lcm_value = sympy.lcm(denominators)

        # Multiply through by denominators to eliminate them
        equation_no_frac = equation_str * lcm_value

        # Rearrange terms
        equation_no_frac = sympy.collect(equation_no_frac, (x, y))
        # Simplify the equation and return the string form
        formula = str(sympy.simplify(equation_no_frac))

        return formula

    def __repr__(self):
        return f"Original : {self.original_input}"

    def __str__(self):
        """
        Return the equation in user-friendly syntax
        """
        equation = self.to_python_syntax()
        return equation.replace('**', '^').replace('*', '')

    def to_python_syntax(self):
        """
        Return the equation in Python syntax
        Probably a waste of time as can use self.equation
        """
        equation = ''
        for powers, coeff in sorted(self.coeff.items(), reverse=True):
            term = ''
            if powers == (2, 0):  # term in x^2
                term += f'{coeff}*x**2'
            elif powers == (1, 1):  # term in xy
                term += f'{coeff}*x*y'
            elif powers == (0, 2):  # term in y^2
                term += f'{coeff}*y**2'
            elif powers == (1, 0):  # term in x
                term += f'{coeff}*x'
            elif powers == (0, 1):  # term in y
                term += f'{coeff}*y'
            elif powers == (0, 0):  # constant term
                term += str(coeff)
            equation += ' + ' + term if coeff > 0 else ' - ' + term[1:]
        return equation.lstrip(' + ')

    def classify(self):
        """
        Classify the conic section based on its coefficients.
        """
        A = self.coefficients.get((2, 0), 0)
        B = self.coefficients.get((1, 1), 0)
        C = self.coefficients.get((0, 2), 0)
        D = self.coefficients.get((1, 0), 0)
        E = self.coefficients.get((0, 1), 0)
        F = self.coefficients.get((0, 0), 0)

        discriminant = B**2 - 4*A*C

        if discriminant < 0:
            return "Circle" if A == C and B == 0 else "Ellipse"
        elif discriminant == 0:
            return "Parabola"
        elif discriminant > 0:
            return "Hyperbola"
        else:
            return "Unknown"

    def get_standard_coefficients(self):
        A = self.coeff.get((2, 0), 0)
        B = self.coeff.get((1, 1), 0) * 2  # This is because the coefficient of xy in the equation is B/2
        C = self.coeff.get((0, 2), 0)
        D = self.coeff.get((1, 0), 0) * 2  # This is because the coefficient of x in the equation is D/2
        E = self.coeff.get((0, 1), 0) * 2  # This is because the coefficient of y in the equation is E/2
        F = self.coeff.get((0, 0), 0)
        return A, B, C, D, E, F

    def plot(self):
        A, B, C, D, E, F = self.get_standard_coefficients()

        x = np.linspace(-10, 10, 400)
        y = np.linspace(-10, 10, 400)
        x, y = np.meshgrid(x, y)

        plt.contour(x, y, (A * x ** 2 + B * x * y + C * y ** 2 + D * x + E * y + F), [1], colors='r')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.title(f"${self.__str__()}$")
        plt.show()

    def compute_transformation(self):
        # Compute the transformation matrix
        pass

    def to_standard_form(self):
        # Apply the transformation matrix to the equation and convert it to standard form
        pass


class Parabola(Conic):
    def __init__(self, equation, coefficients):
        super().__init__(equation)
        # Compute any properties unique to parabolas here

    def plot(self):
        pass  # Implement drawing for parabola


class Ellipse(Conic):
    def __init__(self, equation, coefficients):
        super().__init__(equation)
        # Compute any properties unique to ellipses here

    def plot(self):
        pass  # Implement plotting for ellipse


class Circle(Conic):
    def __init__(self, equation, coefficients):
        super().__init__(equation)
        # Compute any properties unique to circles here

    def plot(self):
        pass  # Implement plotting for circle


class Hyperbola(Conic):
    def __init__(self, equation, coefficients):
        super().__init__(equation)
        # Compute any properties unique to hyperbolas here

    def plot(self):
        pass  # Implement plotting for hyperbola



