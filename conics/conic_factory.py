import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import sympy
from fractions import Fraction
from IPython.display import Math, display
from conics.poly_dictionary import decompose_polynomial
from conics.plotting import *


class Conic:
    def __init__(self, equation: str):
        """
        A class representing a general conic section.
        :param: an equation of a conic section as a string
        """
        self.original_input = equation
        self.equation_str, self.equation = self._to_general_form_conic(equation)
        self.coefficients = decompose_polynomial(self.equation_str)
        self._validate_conic()
        self.A = self.coeff.get((2, 0), 0)
        self.B = self.coeff.get((1, 1), 0)
        self.C = self.coeff.get((0, 2), 0)
        self.D = self.coeff.get((1, 0), 0)
        self.E = self.coeff.get((0, 1), 0)
        self.F = self.coeff.get((0, 0), 0)
        self.coeff_matrix = self.create_coeff_matrix()
        self.quad_matrix = self.create_quad_matrix()
        self.type = self.classify()
        self.expression = None
        self.history = []

    @classmethod
    def create(cls, equation: str):
        """
        Creates an instance of a Conic subclass based on the type of conic section.
        -----
        Raises;
        ValueError if the type of the conic section cannot be determined.
        """
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
            raise ValueError('Cannot instantiate this type of the conic section')

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
        equation = str(equation).replace('**', '^').replace('*', '').replace('sqrt', '√')
        equation += ' = 0'
        return equation

    def _to_general_form_conic(self, equation: str):
        """
        Takes the string, works out the LCM, collects like terms
        then reduces the resultant as far as possible.
        :return: General form equation of a conic as string and sympy expression
        """
        x, y = sympy.symbols('x y')
        equation_str = sympy.sympify(equation, locals={"sqrt": sympy.sqrt})

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
        formula = sympy.expand(formula)     # General form requires ZERO arbitrary factorisation
        formula_str = str(formula)

        return formula_str, formula

    def _validate_conic(self):
        if not any(sum(key) == 2 for key in self.coefficients) \
                or any(sum(key) > 2 for key in self.coefficients):
            raise ValueError("The polynomial does not represent a conic section.")

    def create_coeff_matrix(self):
        """
        Creates the matrix of coefficients,
        """
        matrix = np.array([[self.A, self.B / 2, self.D / 2],
                           [self.B / 2, self.C, self.E / 2],
                           [self.D / 2, self.E / 2, self.F]])
        return matrix

    def create_quad_matrix(self):
        """
        Creates a 2x2 sub-matrix from the coefficient matrix
        """
        matrix = np.array([[self.A, self.B / 2],
                           [self.B / 2, self.C]])
        return matrix

    def print_matrices(self):
        print(f"Matrix of coefficients:\n{self.coeff_matrix}\n\nQuadratic Matrix:\n{self.quad_matrix}")

    def classify(self):
        """
        Classify the conic section based on its coefficients.
        Uses the determinant of the 2x2 sub-matrix (delta) and the determinant of the 3x3 matrix (Delta).
        """
        # Calculate the determinant of the 2x2 quadratic-matrix
        delta = np.linalg.det(self.quad_matrix)

        # TODO: Methods remain, useful for classification however, unused
        # Calculate trace of the 2x2 sub-matrix
        tau = np.trace(self.quad_matrix)

        # Calculate the determinant of the 3x3 matrix. - DEGENERATE CONICS?
        Delta = np.linalg.det(self.coeff_matrix)

        if Delta == 0:
            return "Degenerate"
        elif delta == 0:
            return "Parabola"
        elif delta < 0:
            return "Hyperbola"
        elif delta > 0:
            if self.A == self.C and self.B == 0:
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
        display(self.expression.as_expr())

    def save_as_sympy(self, rational=False, return_expr=False):
        """
        Saves the polynomial as a sympy expression and returns it.
        """
        self._convert_to_sympy_expression(rational, force=True)
        if return_expr:
            return self.expression.as_expr()  # Return as a sympy expression
        else:
            return self.expression  # Return as a sympy.Poly object

    def draw(self):
        """
        Method to quickly plot an instance of a conic.
        Superseded by more advanced plotting operations in subclasses
        """
        x = np.linspace(-10, 10, 400)
        y = np.linspace(-10, 10, 400)
        x, y = np.meshgrid(x, y)

        plt.contour(x, y, (self.A * x ** 2 + self.B * x * y + self.C * y ** 2
                           + self.D * x + self.E * y + self.F), [1], colors='r')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.axhline(0, color='gray', linewidth=0.5)
        plt.axvline(0, color='gray', linewidth=0.5)
        plt.title(f"${self.__str__()}$")
        plt.show()


class Parabola(Conic):
    """
    A class that represents a Parabola, extending the Conic class.
    """
    def __init__(self, equation):
        super().__init__(equation)
        self.history = []
        self.standard_form = False

    @property
    def orientation(self):
        return self.get_orientation()

    @property
    def axis(self):
        return self.compute_axis()

    @property
    def vertex(self):
        return self.compute_vertex()

    @property
    def a(self):
        if not self.standard_form:
            raise ValueError("The parabola is not in standard form")
        return - self.D / (4 * self.C)

    @property
    def focus(self):
        if not self.standard_form:
            raise ValueError("The parabola is not in standard form")
        return self.a, 0

    @property
    def directrix(self):
        if not self.standard_form:
            raise ValueError("The parabola is not in standard form")
        return - self.a

    @property
    def latus_rectum(self):
        if not self.standard_form:
            raise ValueError("The parabola is not in standard form")
        return 4 * self.a

    def get_info(self):
        print(f"{self.__repr__()}\nType: {self.type}\nCoefficients: {self.coeff}"
              f"\nGeneral Form: {self}\n")
        self.print_matrices()
        print(f"\nOrientation: {self.get_orientation()}")
        print(f"Axis of symmetry: {str(self.axis)}")
        self.plot_parabola()

    def get_orientation(self):
        """
        Determine the orientation of the parabola in E^2.

        The orientation is computed based on the coefficients of the parabola equation.
        If the parabola is rotated, the method also returns the angle of rotation in radians.

        Returns
        -------
        str, str or float
            If parabola is vertical or horizontal, returns a tuple of two strings:
            the first string represents the orientation ('vertical' or 'horizontal'),
            the second string indicates the direction ('positive' or 'negative').

            If parabola is rotated, returns a tuple where the first element is 'rotated',
            and the second element is the angle of rotation in radians.
        """
        if self.A != 0 and self.C == 0:
            if self.A > 0:
                return "vertical", "positive"
            else:
                return "vertical", "negative"
        elif self.A == 0 and self.C != 0:
            if self.C > 0:
                if self.D > 0:
                    return "horizontal", "negative"  # modification here
                else:  # when self.D =< 0
                    return "horizontal", "positive"  # and here
            else:
                if self.D > 0:
                    return "horizontal", "positive"
                else:  # when self.D =< 0
                    return "horizontal", "negative"
        elif self.A != 0 and self.C != 0:
            theta = 0.5 * sympy.atan2(self.B, (self.A - self.C))
            return "rotated", theta

    def compute_axis(self):
        """
        Computes the axis of symmetry

        For a non-rotated parabola, the axis is computed using the coefficients of the equation.
        For a rotated parabola, the axis is computed using the derivatives of the general equation.

        Returns
        -------
        The equation of the axis in the form of a sympy Rational number for non-rotated parabolas.
        For rotated parabolas, returns a sympy Add object representing the equation of the axis.
        """
        if not self.orientation[0] == 'rotated':
            if self.orientation[0] == 'vertical':
                axis = sympy.Rational(-self.D / (2 * self.A)).limit_denominator(1000000)
            else:
                axis = sympy.Rational(-self.E / (2 * self.C)).limit_denominator(100000)
        else:
            x, y = sympy.symbols('x y')
            gen_eqn = self.save_as_sympy(return_expr=True)

            # Compute the derivatives
            gen_x = sympy.diff(gen_eqn, x)
            gen_y = sympy.diff(gen_eqn, y)

            # Factorize the quadratic part
            alpha = sympy.sqrt(self.A)
            if self.B >= 0:
                beta = sympy.sqrt(self.C)
            else:
                beta = -sympy.sqrt(self.C)
            axis = alpha * gen_x + beta * gen_y

        return axis

    def compute_vertex(self):
        """
        Computes the vertex.

        For a non-rotated parabola, the vertex is computed using the coefficients of the equation.
        The method will handle both vertical and horizontal parabolas.
        For a rotated parabola, a different method is employed, and None is returned for both h and k.
        """
        if self.orientation[0] != 'rotated':  # handle vertical and horizontal
            if self.orientation[0] == 'vertical':  # Parabola opens up or down
                a = self.A / -self.E
                d = self.D / -self.E
                f = self.F / -self.E
                h = self._rational_or_radical(-d / (2 * a))
                k = self._rational_or_radical(a * h ** 2 + d * h + f)
            else:  # Parabola opens to the right or left
                c = self.C / -self.D
                e = self.E / -self.D
                f = self.F / -self.D
                k = self._rational_or_radical(-e / (2 * c))
                h = self._rational_or_radical(c * k ** 2 + e * k + f)
        else:  # handle rotated case
            h = None
            k = None

        return h, k

    def _rational_or_radical(self, x):
        """
        Convert x to a rational number if it is a float, otherwise leave it as it is.
        """
        if isinstance(x, (float, int, sympy.Float)):
            return sympy.Rational(x).limit_denominator(100000)
        else:
            return x

    def rotate_parabola(self, rational=False, display=False):
        """
        Rotate the parabola until orientation, "horizontal, positive" or "standard form".

        This method uses a rotation matrix to rotate the parabola. The rotation angle is
        determined based on the parabola's orientation.

        This function is recursive: if the resulting parabola is not a "horizontal, positive"
        one after rotation, this function will be called again until the desired orientation
        is achieved.

        The rotation operation is recorded in the object's history.
        """
        self.history.append(f"Equation: {str(self)}")
        orientation, rotation_angle = self.get_orientation()

        if orientation == "vertical":
            if rotation_angle == "positive":
                rotation_angle = sympy.pi * 3 / 2  # 270 degrees
            elif rotation_angle == "negative":
                rotation_angle = sympy.pi / 2  # 90 degrees
        elif orientation == "horizontal":
            rotation_angle = sympy.pi if rotation_angle == "negative" else 0

        if orientation == "rotated":
            epsilon = 1e-10
            # We use a new variable here to avoid the issue
            rotated_angle = rotation_angle
            if abs(abs(rotated_angle) - sympy.pi / 4) > epsilon:  # Not close to 45 degrees
                rotated_angle = sympy.pi / 2 - rotated_angle
            rotation_angle = rotated_angle

        # Create Rotation Matrix
        R = np.array([[sympy.cos(rotation_angle), sympy.sin(rotation_angle), 0],
                      [-sympy.sin(rotation_angle), sympy.cos(rotation_angle), 0],
                      [0, 0, 1]])

        # Rotate parabola
        self.coeff_matrix = np.transpose(R) @ self.coeff_matrix @ R
        threshold = 1e-14
        self.coeff_matrix = np.where(abs(self.coeff_matrix) < threshold, 0, self.coeff_matrix)

        if rational:
            for i in range(self.coeff_matrix.shape[0]):
                for j in range(self.coeff_matrix.shape[1]):
                    self.coeff_matrix[i, j] = sympy.Rational(self.coeff_matrix[i, j])
                # Now all entries in coeff_matrix are SymPy Rationals

        # Update state from the new matrix and record the new state in history
        self.history.append(f"Performed rotation by {sympy.N(rotation_angle * 180/sympy.pi, 4)} degrees CCW")
        self._update_from_matrix()
        self.record_state()

        if display:
            self.symbolic_rotation()

        # After rotating, if the parabola isn't oriented as 'horizontal, positive',
        # recursively rotate again
        orientation, rotation_angle = self.get_orientation()
        if orientation == "vertical" or (orientation == "horizontal" and rotation_angle == "negative"):
            self.rotate_parabola(rational)

    def symbolic_rotation(self):
        # Define your rotation matrix R
        theta = sympy.symbols('theta')
        R = sympy.Matrix([
            [sympy.cos(theta), sympy.sin(theta), 0],
            [-sympy.sin(theta), sympy.cos(theta), 0],
            [0, 0, 1]])

        # Suppose your coeff_matrix is C (a 3x3 matrix)
        A, B, C, D, E, F = sympy.symbols('A B C D E F')
        M = sympy.Matrix([
            [A, B / 2, D / 2],
            [B / 2, C, E / 2],
            [D / 2, E / 2, F]])

        # Display the symbolic operation
        display(Math('\\textbf{R}^T \cdot \\textbf{M} \cdot \\textbf{R}'))
        display(Math(f'{sympy.latex(R.transpose())} \cdot {sympy.latex(M)} \cdot {sympy.latex(R)}'))

    def translate_origin(self, rational=False, display=False):
        """
        Translate the parabola so vertex is at the origin.

        This method uses an affine transformation.

        The translation operation is recorded in the object's history. After translation, if the
        vertex of the parabola is at the origin and its orientation is "horizontal, positive",
        the parabola is considered to be in its standard form.
        """
        h, k = self.vertex[0], self.vertex[1]
        original = h, k

        # Translation Matrix
        T = np.array([[1, 0, h],
                      [0, 1, k],
                      [0, 0, 1]])

        self.coeff_matrix = T.T @ self.coeff_matrix @ T
        threshold = 1e-14
        self.coeff_matrix = np.where(abs(self.coeff_matrix) < threshold, 0, self.coeff_matrix)

        if rational:
            for i in range(self.coeff_matrix.shape[0]):
                for j in range(self.coeff_matrix.shape[1]):
                    self.coeff_matrix[i, j] = sympy.Rational(self.coeff_matrix[i, j])
                # Now all entries in coeff_matrix are SymPy Rationals

        # Update state from the new matrix and record the new state in history
        self.history.append(f"Affine transformation of vertex {original} to the origin")
        self._update_from_matrix()
        self.record_state()

        if display:
            self.symbolic_translation()

        # Update standard_form flag
        if self.vertex == (0, 0) and self.orientation == ('horizontal', 'positive'):
            self.standard_form = True

    def symbolic_translation(self):
        h, k = sympy.symbols('h k')
        T = sympy.Matrix([[1, 0, h],
                          [0, 1, k],
                          [0, 0, 1]])

        A, B, C, D, E, F = sympy.symbols('A B C D E F')
        M = sympy.Matrix([
            [A, B / 2, D / 2],
            [B / 2, C, E / 2],
            [D / 2, E / 2, F]])

        # Display the symbolic operation
        display(Math('\\textbf{T}^T \cdot \\textbf{M} \cdot \\textbf{T}'))
        display(Math(f'{sympy.latex(T.transpose())} \cdot {sympy.latex(M)} \cdot {sympy.latex(T)}'))

    def _update_coefficients(self):
        self.A = self.coeff_matrix[0, 0]
        self.B = 2 * self.coeff_matrix[0, 1]  # Remember we stored B/2 in the matrix
        self.C = self.coeff_matrix[1, 1]
        self.D = 2 * self.coeff_matrix[0, 2]  # Remember we stored D/2 in the matrix
        self.E = 2 * self.coeff_matrix[1, 2]  # Remember we stored E/2 in the matrix
        self.F = self.coeff_matrix[2, 2]

    def _update_from_matrix(self):
        # Update coefficients from the matrix
        self._update_coefficients()

        # Create the new coefficient dictionary
        self.coefficients = {
            (2, 0): self.A,
            (1, 1): self.B,
            (0, 2): self.C,
            (1, 0): self.D,
            (0, 1): self.E,
            (0, 0): self.F
        }

        # Update the equation string
        self.expression = self.save_as_sympy(return_expr=True)
        self.equation = str(self.expression)

    def record_state(self):
        """
        This method records the current state of the parabola in the history of the object.
        """
        # Store the current state in history
        self.history.append({
            'Coefficients': self.coefficients,
            'Matrix': self.coeff_matrix,
            'Equation': str(self),
            'Orientation': self.orientation,
            'Vertex': self.vertex
        })

    def print_history(self):
        """
        This method prints the history of the transformations applied to the parabola
        """
        print(f"Original input:\n{self.original_input}\n")
        for item in self.history:
            if isinstance(item, str):
                print(f"----------\n{item}\n")
            elif isinstance(item, dict):
                for key, value in item.items():
                    print(f"{key} :\n{value}\n")

    def plot_parabola(self, x_range=None, y_range=None):
        plot_parabola(self, x_range, y_range)

    def plot_standard(self, x_range=None, y_range=None):
        parabola_standard(self, x_range, y_range)


class Circle(Conic):
    def __init__(self, equation):
        super().__init__(equation)
        self.history = []

    @property
    def radius(self):
        return self.compute_radius()

    @property
    def centre(self):
        return self.compute_centre()

    @property
    def area(self):
        return self.compute_area()

    @property
    def circumference(self):
        return self.compute_circumference()

    def get_info(self, x_range=None, y_range=None):
        print(f"{self.__repr__()}\nType: {self.type}\nCoefficients: {self.coeff}"
              f"\nGeneral Form: {self}\n\nRadius: {self.radius:.2f}\n\n"
              f"Area: {self.area[0]}\napprox {self.area[1]} sq. units"
              f"\n\nCircumference: {self.circumference[0]}\napprox {self.circumference[1]} units\n")
        self.print_matrices()
        self.plot_circle(x_range, y_range)

    def compute_centre(self):
        """Compute and return the center of the circle."""
        return -self.D / 2, -self.E / 2

    def compute_radius(self):
        """Compute and return the radius of the circle."""
        return sympy.sqrt((self.D ** 2 / 4) + (self.E ** 2 / 4) - self.F)

    def compute_area(self):
        return sympy.pi * self.radius ** 2, sympy.N(sympy.pi * self.radius ** 2, 6)

    def compute_circumference(self):
        return 2 * sympy.pi * self.radius, sympy.N(2 * sympy.pi * self.radius, 6)

    def plot_circle(self, x_range=None, y_range=None):
        plot_circle(self, x_range, y_range)


class Ellipse(Conic):
    def __init__(self, equation):
        super().__init__(equation)
        # Compute any properties unique to ellipses here

    def draw(self):
        pass  # Implement plotting for ellipse


class Hyperbola(Conic):
    def __init__(self, equation):
        super().__init__(equation)
        # Compute any properties unique to hyperbolas here

    def draw(self):
        pass  # Implement plotting for hyperbola



