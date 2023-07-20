import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import sympy
from abc import abstractmethod
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
        elif conic_type == "Hyperbola" or conic_type == "Rectangular Hyperbola":
            return Hyperbola(equation)
        else:
            raise ValueError('Cannot instantiate this type of the conic section')

    @abstractmethod
    def _get_rotation_angle(self):
        """Calculate and return the rotation angle for the conic section."""
        pass

    @abstractmethod
    def _get_translation_point(self):
        """Translation point will be the vertex for Parabolas and centre for ellipse/hyperbolas"""
        pass

    @abstractmethod
    def _standard_form_flag(self):
        """Updates the standard form flag in the parent class"""
        pass

    @property
    def coeff(self):
        return self.coefficients

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
        matrix = sympy.Matrix([[self.A, sympy.Rational(self.B, 2), sympy.Rational(self.D, 2)],
                               [sympy.Rational(self.B, 2), self.C, sympy.Rational(self.E, 2)],
                               [sympy.Rational(self.D, 2), sympy.Rational(self.E, 2), self.F]])
        return matrix

    def create_quad_matrix(self):
        """
        Creates a 2x2 sub-matrix from the coefficient matrix
        """
        matrix = sympy.Matrix([[self.A, sympy.Rational(self.B, 2)],
                               [sympy.Rational(self.B, 2), self.C]])
        return matrix

    def print_matrices(self):
        print("Matrix of coefficients:")
        display(self.coeff_matrix)
        print("\nQuadratic Matrix:")
        display(self.quad_matrix)

    def classify(self):
        """
        Classify the conic section based on its coefficients.
        Uses the determinant of the 2x2 sub-matrix (delta) and the determinant of the 3x3 matrix (Delta).
        """
        # Calculate the determinant of the 2x2 quadratic-matrix
        delta = self.quad_matrix.det()

        # Calculate trace of the 2x2 sub-matrix
        tau = self.quad_matrix.trace()

        # Calculate the determinant of the 3x3 matrix. - DEGENERATE CONICS?
        Delta = self.coeff_matrix.det()

        if Delta == 0:
            return "Degenerate"
        elif delta == 0:
            return "Parabola"
        elif delta < 0:
            if tau == 0:
                return "Rectangular Hyperbola"
            else:
                return "Hyperbola"
        elif delta > 0:
            if self.A == self.C and self.B == 0:
                return "Circle"
            else:
                return "Ellipse"
        else:
            return "Unknown"

    def _clean_coeff_matrix(self, threshold, tolerance, rational):
        for i in range(3):
            for j in range(3):
                if abs(self.coeff_matrix[i, j]) < threshold:
                    self.coeff_matrix[i, j] = 0
                if rational:
                    self.coeff_matrix[i, j] = sympy.nsimplify(self.coeff_matrix[i, j], tolerance=tolerance)

    def _update_coefficients(self):
        self.A = self.coeff_matrix[0, 0]
        self.B = 2 * self.coeff_matrix[0, 1]
        self.C = self.coeff_matrix[1, 1]
        self.D = 2 * self.coeff_matrix[0, 2]
        self.E = 2 * self.coeff_matrix[1, 2]
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

    def rotate(self, display=False, rational=True, threshold=1e-14, tolerance=0.001):
        """
        This method uses a rotation matrix to rotate the conic section.
        The rotation angle is determined using an abstract method unique within each subclass

        The rotation operation is recorded in the object's history.
        """
        self.history.append(f"Equation: {str(self)}")
        if self.type == 'Parabola':
            self.orientation_history.append(self.coefficients)

        rotation_angle = self._get_rotation_angle()

        R = sympy.Matrix([[sympy.cos(rotation_angle), sympy.sin(rotation_angle), 0],
                          [-sympy.sin(rotation_angle), sympy.cos(rotation_angle), 0],
                          [0, 0, 1]])

        self.coeff_matrix = R.transpose() * self.coeff_matrix * R

        # if the absolute value of an element in coeff_matrix is less than threshold, we set it to 0
        self._clean_coeff_matrix(threshold, tolerance, rational)

        self.history.append(f"Performed rotation by {round(float(sympy.deg(rotation_angle)), 4)} degrees CCW")
        self._update_from_matrix()
        self._standard_form_flag()

        # Update History
        self.record_state()

        if display:
            self._symbolic_rotation()

    def translate(self, rational=False, display=False, threshold=1e-14, tolerance=0.001):
        """
        Translates the vertex/centre of the conic to the origin.

        The translation operation is recorded in the object's history.
        """
        h, k, original = self._get_translation_point()

        # Translation Matrix
        T = sympy.Matrix([[1, 0, h],
                          [0, 1, k],
                          [0, 0, 1]])

        self.coeff_matrix = T.T * self.coeff_matrix * T

        # if the absolute value of an element in coeff_matrix is less than threshold, we set it to 0
        self._clean_coeff_matrix(threshold, tolerance, rational)

        # Update state from the new matrix and record the new state in history
        point = 'vertex' if self.type == "Parabola" else 'centre'
        self.history.append(f"Transformation of {point} {original} to the origin")
        self._update_from_matrix()
        self._standard_form_flag()

        # Update history
        self.record_state()

        if display:
            self._symbolic_translation()

    def _symbolic_rotation(self):
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
        print("We form a rotation matrix R:")
        display(Math('\\textbf{R}^T \cdot \\textbf{M} \cdot \\textbf{R}'))
        print("Such that,")
        display(Math(f'{sympy.latex(R.transpose())} \cdot {sympy.latex(M)} \cdot {sympy.latex(R)}'))

    def _symbolic_translation(self):
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
        point = 'vertex' if self.type == "Parabola" else 'centre'
        print(f"We form a translation matrix T where (h, k) are the coordinates (x, y) of the {point}")
        display(Math('\\textbf{T}^T \cdot \\textbf{M} \cdot \\textbf{T}'))
        print("Such that")
        display(Math(f'{sympy.latex(T.transpose())} \cdot {sympy.latex(M)} \cdot {sympy.latex(T)}'))

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

    def save_as_sympy(self, rational=False, return_expr=False):
        """
        Saves the polynomial as a sympy expression and returns it.
        """
        self._convert_to_sympy_expression(rational, force=True)
        if return_expr:
            return self.expression.as_expr()  # Return as a sympy expression
        else:
            return self.expression  # Return as a sympy.Poly object

    def display(self):
        """
        Displays the sympy expression of the polynomial using IPython's display system.
        """
        self._convert_to_sympy_expression()
        display(self.expression.as_expr())

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
        self.orientation_history = []
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
        if self.orientation[0] == 'Rotated':
            angle_in_degrees = sympy.N(sympy.deg(self.orientation[1]), 5)
            print(f"Orientation: {self.orientation[0], angle_in_degrees} degrees\nradians: {self.orientation[1]}\n")
        else:
            print(f"Orientation: {self.orientation}\n")
        print(f"Axis of symmetry:")
        display(self.axis)
        self.plot()

    def get_orientation(self):
        """
        Determine the orientation of the parabola in E^2.

        The orientation is computed based on the coefficients of the parabola equation.
        If the parabola is rotated, the method also returns the angle of rotation in radians.
        """
        if self.A != 0 and self.C == 0:
            if self.A > 0:
                return "Vertical", "Positive"
            else:
                return "Vertical", "Negative"

        elif self.A == 0 and self.C != 0:
            if self.C > 0:
                if self.D > 0:
                    return "Horizontal", "Negative"  # modification here
                else:  # when self.D =< 0
                    return "Horizontal", "Positive"  # and here
            else:
                if self.D > 0:
                    return "Horizontal", "Positive"
                else:  # when self.D =< 0
                    return "Horizontal", "Negative"

        elif self.A != 0 and self.C != 0:
            theta = sympy.Rational(1, 2) * sympy.atan2(self.B, (self.A - self.C))

            return "Rotated", float(theta)

    def compute_axis(self):
        """
        Computes the axis of symmetry.

        For a non-rotated parabola, the axis is computed using the coefficients of the equation.
        For a rotated parabola, the axis is computed using the derivatives of the general equation.

        Returns
        -------
        The equation of the axis in the form of a sympy Equality object for all parabolas.
        """
        x, y = sympy.symbols('x y')
        if not self.orientation[0] == 'Rotated':
            if self.orientation[0] == 'Vertical':
                # The axis of symmetry for a vertical parabola is a vertical line (x = h)
                h = sympy.Rational(-self.D / (2 * self.A)).limit_denominator(100000)
                axis = sympy.Eq(x, h)
            else:
                # The axis of symmetry for a horizontal parabola is a horizontal line (y = k)
                k = sympy.Rational(-self.E / (2 * self.C)).limit_denominator(100000)
                axis = sympy.Eq(y, k)
        else:
            x, y = sympy.symbols('x y')

            # Save the general equation in the symbolic form
            gen_eqn = self.save_as_sympy(return_expr=True)

            # Compute the derivatives
            gen_x = sympy.diff(gen_eqn, x)
            gen_y = sympy.diff(gen_eqn, y)

            # Factorize the quadratic part
            alpha = sympy.sqrt(self.A)
            beta = sympy.sqrt(self.C) if self.B >= 0 else -sympy.sqrt(self.C)

            # Form the equation alpha * gen_x + beta * gen_y = 0
            equation = alpha * gen_x + beta * gen_y

            # Rearrange the equation in the form y = mx + b
            equation_y = sympy.solve(equation, y)[0]

            # Simplify the equation
            simplified_equation_y = sympy.simplify(equation_y)

            # Form an equation object
            axis = sympy.Eq(y, simplified_equation_y)

        return axis

    def compute_vertex(self):
        """
        Computes the vertex.

        For a non-rotated parabola, the vertex is computed using the coefficients of the equation.
        The method will handle both vertical and horizontal parabolas.
        For a rotated parabola, a different method is employed, and None is returned for both h and k.
        """
        if self.orientation[0] != 'Rotated':  # handle vertical and horizontal
            if self.orientation[0] == 'Vertical':  # Parabola opens up or down
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

    def get_original_vertex(self, to_float=False, tolerance=0.001):
        """
        This method only applies to once rotated parabolas
        :param: default rational number True
        :return:
        """
        # Fetch the vertex of the rotated parabola
        vertex_rotated = self.vertex

        # Extract the original coefficients from the history (state before rotation)
        original_coefficients = self.orientation_history[0]

        # Use a dictionary comprehension to extract the original A, B, and C values
        original_A = original_coefficients.get((2, 0), 0)
        original_B = original_coefficients.get((1, 1), 0)
        original_C = original_coefficients.get((0, 2), 0)

        # Calculate theta using the original A, B, and C values
        theta = sympy.Rational(1, 2) * sympy.atan2(original_B, (original_A - original_C))
        # Correct for the rotation by subtracting pi/2
        theta -= sympy.pi / 2
        theta = sympy.nsimplify(theta, tolerance=tolerance)

        # Unpack the rotated vertex coordinates
        x_rotated, y_rotated = vertex_rotated

        # Rotate the vertex coordinates back to the original coordinate system
        x_vertex_expr = x_rotated * sympy.cos(theta) - y_rotated * sympy.sin(theta)
        y_vertex_expr = x_rotated * sympy.sin(theta) + y_rotated * sympy.cos(theta)

        # Calculate the rational representation with limited denominator size
        x_vertex = sympy.nsimplify(x_vertex_expr, tolerance=tolerance)
        y_vertex = sympy.nsimplify(y_vertex_expr, tolerance=tolerance)

        if to_float:
            x_vertex = float(x_vertex)
            y_vertex = float(y_vertex)

        return x_vertex, y_vertex

    def _get_rotation_angle(self):
        orientation, rotation_angle = self.get_orientation()

        if orientation == "Vertical":
            if rotation_angle == "Positive":
                rotation_angle = sympy.pi * 3 / 2  # 270 degrees
            elif rotation_angle == "Negative":
                rotation_angle = sympy.pi / 2  # 90 degrees
        elif orientation == "Horizontal":
            rotation_angle = sympy.pi if rotation_angle == "Negative" else 0

        if orientation == "Rotated":
            epsilon = 1e-10
            rotated_angle = rotation_angle
            if abs(abs(rotated_angle) - sympy.pi / 4) > epsilon:  # Not close to 45 degrees
                rotated_angle = sympy.pi / 2 - rotated_angle
            rotation_angle = rotated_angle

        return rotation_angle

    def _get_translation_point(self):
        h, k = self.vertex[0], self.vertex[1]
        original = h, k
        return h, k, original

    def _standard_form_flag(self):
        # Update standard_form flag
        if self.vertex == (0, 0) and self.orientation == ('Horizontal', 'Positive'):
            self.standard_form = True

    def _rational_or_radical(self, x):
        """
        Convert x to a rational number if it is a float, otherwise leave it as it is.
        """
        if isinstance(x, (float, int, sympy.Float)):
            return sympy.Rational(x).limit_denominator(100000)
        else:
            return x

    def record_state(self):
        """
        This method records the current state of the parabola in the history of the object.
        """
        self.history.append({
            'Coefficients': self.coefficients,
            'Matrix': self.coeff_matrix,
            'Equation': str(self),
            'Orientation': self.orientation,
            'Vertex': self.vertex,
            'Standard Form?': self.standard_form
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

    def plot(self, x_range=None, y_range=None):
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
        self.plot(x_range, y_range)

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

    def rotate(self, *args, **kwargs):
        print("Rotation operation is not applicable for circles.")

    def _get_translation_point(self):
        """Abstract method of translate"""
        h, k = self.centre[0], self.centre[1]
        original = h, k
        return h, k, original

    def plot(self, x_range=None, y_range=None):
        plot_circle(self, x_range, y_range)


class Ellipse(Conic):
    def __init__(self, equation):
        super().__init__(equation)
        self.standard_form = False

    @property
    def centre(self):
        return self.get_centre()

    @property
    def orientation(self):
        """
        Return the orientation of the ellipse: 'Horizontal', 'Vertical' or 'Rotated'.
        """
        # Calculate the angle in radians
        theta = sympy.Rational(1, 2) * sympy.atan(sympy.Rational(self.B, (self.A - self.C)))

        # Use sympy's equivalence checking for comparison
        if self.B != 0:
            return 'Rotated', theta
        else:
            if self.A > self.C:
                return 'Vertical', 0
            elif self.A < self.C:
                return 'Horizontal', 0

    @property
    def semimajor_axis(self):
        if not self.standard_form:
            raise ValueError("Semi-major axis can only be computed for ellipses in standard form")
        a = self.A
        return sympy.sqrt(1 / a)

    @property
    def semiminor_axis(self):
        if not self.standard_form:
            raise ValueError("Semi-minor axis can only be computed for ellipses in standard form")
        b = self.C
        return sympy.sqrt(1 / b)

    @property
    def eccentricity(self):
        if not self.standard_form:
            raise ValueError("Eccentricity can only be computed for ellipses in standard form")
        # eccentricity e = sqrt(1 - (b/a)^2), a > b
        a, b = self.semimajor_axis, self.semiminor_axis
        e = sympy.sqrt(1 - (b/a)**2)
        return e

    @property
    def vertices(self):
        if not self.standard_form:
            raise ValueError("Vertices can only be computed for ellipses in standard form")
        # Vertices are on the semi-major axis, a distance of 'a' from the center
        a, b = self.semimajor_axis, self.semiminor_axis
        major_vertices = [(self.centre[0] - a, self.centre[1]), (self.centre[0] + a, self.centre[1])]
        minor_vertices = [(self.centre[0], self.centre[1] - b), (self.centre[0], self.centre[1] + b)]
        return major_vertices, minor_vertices

    @property
    def foci(self):
        if not self.standard_form:
            raise ValueError("Foci can only be computed for ellipses in standard form")
        # Foci are on the semi-major axis, a distance of +/- 'ae' from the center
        a, e = self.semimajor_axis, self.eccentricity
        return [(self.centre[0] - a*e, self.centre[1]), (self.centre[0] + a*e, self.centre[1])]

    @property
    def directrices(self):
        if not self.standard_form:
            raise ValueError("Directrices can only be computed for ellipses in standard form")
        # Directrices are lines parallel to the minor axis, a distance of +/- 'a/e' from the center
        a, e = self.semimajor_axis, self.eccentricity
        x = sympy.symbols('x')
        return [sympy.Eq(x, (self.centre[0] - a / e)), sympy.Eq(x, (self.centre[0] + a / e))]

    def get_info(self):
        print(f"{self.__repr__()}\nType: {self.type}\nCoefficients: {self.coeff}"
              f"\nGeneral Form: {self}\n")
        self.print_matrices()
        if self.orientation[0] == 'Rotated':
            angle_in_degrees = sympy.N(sympy.deg(self.orientation[1]), 5)
            print(f"Orientation: {self.orientation[0], angle_in_degrees} degrees\nradians: {self.orientation[1]}\n")
        else:
            print(f"Orientation: {self.orientation[0]}\n")
        print(f"Centre: {self.centre}")
        print(f"Equation of line through centre & semi-major axis")
        display(self.semimajor_axis_line())
        self.plot()

    def semimajor_axis_line(self):
        """
        Return the equation of the line passing through the semi-major axis of the ellipse.
        """
        x, y = sympy.symbols('x y')

        # Calculate the slope from the rotation angle
        if self.orientation[0] == 'Rotated':
            m = sympy.tan(self.orientation[1])
            # Calculate the y-intercept from the y-coordinate of the center
            c = self.centre[1] - m * self.centre[0]
            eqn = sympy.nsimplify(sympy.Eq(y, m * x + c), 0.0001)
            simplified_eq = sympy.Eq(eqn.lhs, sympy.trigsimp(eqn.rhs))
            return simplified_eq

        elif self.orientation[0] == 'Horizontal':
            return sympy.Eq(y, self.centre[1])

        else:  # vertical
            return sympy.Eq(x, self.centre[0])

    def get_centre(self):
        """
        The center of the ellipse is the point where the derivatives w.r.t x and y are zero.
        It's a stationary point of the ellipse
        :return: tuple(x_coordinate, y_coordinate)
        """
        # Define the symbols
        x, y = sympy.symbols('x y')

        # Define the system of equations
        eq1 = sympy.Eq(2 * self.A * x + self.B * y + self.D, 0)
        eq2 = sympy.Eq(self.B * x + 2 * self.C * y + self.E, 0)

        # Solve the system
        solution = sympy.solve((eq1, eq2), (x, y))

        center_x = solution[x]
        center_y = solution[y]

        return center_x, center_y

    def _get_rotation_angle(self):
        """Abstract method of rotate"""
        orientation, rotation_angle = self.orientation[0], self.orientation[1]

        if orientation == "Vertical":
            rotation_angle = sympy.pi / 2  # 90

        if orientation == "Rotated":
            if rotation_angle < 0:
                rotation_angle = (abs(rotation_angle)).evalf()
            else:
                rotation_angle = (sympy.pi - rotation_angle).evalf()

        return rotation_angle

    def _get_translation_point(self):
        """Abstract method of translate"""
        h, k = self.centre[0], self.centre[1]
        original = h, k
        return h, k, original

    def _standard_form_flag(self):
        """Abstract method of rotate and translate"""
        # Update standard_form flag
        if self.centre == (0, 0) and self.orientation[0] == 'horizontal':
            self.standard_form = True
            # Normalise the constant term in coefficient matrix
            self.coeff_matrix = abs(1/self.F) * self.coeff_matrix
            # Update coefficients
            self._update_from_matrix()

    def record_state(self):
        """
        This method records the current state of the parabola in the history of the object.
        """
        self.history.append({
            'Coefficients': self.coefficients,
            'Matrix': self.coeff_matrix,
            'Equation': str(self),
            'Orientation': self.orientation,
            'Centre': self.centre,
            'Standard Form?': self.standard_form
            })

    def print_history(self):
        """
        This method prints the history of the transformations applied to the ellipse
        """
        print(f"Original input:\n{self.original_input}\n")
        for item in self.history:
            if isinstance(item, str):
                print(f"----------\n{item}\n")
            elif isinstance(item, dict):
                for key, value in item.items():
                    print(f"{key} :\n{value}\n")

    def plot(self, x_range=None, y_range=None):
        plot_ellipse(self, x_range, y_range)

    def plot_standard(self, x_range=None, y_range=None, decimal_places=4):
        if self.standard_form:
            plot_standard(self, x_range, y_range)
            sympy.init_printing()

            # Eccentricity
            print("------")
            print("Eccentricity: \n")
            print("Approximation: ", sympy.N(self.eccentricity, decimal_places))
            print("Exact: ", end="")
            display(self.eccentricity)

            # Foci
            print("\n------")
            print("Foci: \n")
            for element in self.foci:
                rounded_element = tuple(sympy.N(val, decimal_places) for val in element)
                print("Approximation: ", rounded_element)
                print("Exact: ", end="")
                display(element)

            # Directrices
            print("\n------")
            print("Directrices: \n")
            for element in self.directrices:
                # Get left and right sides of the equation
                lhs = element.lhs
                rhs = element.rhs
                # Approximate the right side
                rounded_rhs = sympy.N(rhs, decimal_places)
                # Print the approximation
                print(f"Approximation: {lhs} = {rounded_rhs}")
                # Display the exact value
                print("Exact: ", end="")
                display(element)

            # Vertices
            print("\n------")
            print("Vertices: \n")
            for sub_list in self.vertices:
                for element in sub_list:
                    rounded_element = tuple(sympy.N(val, decimal_places) for val in element)
                    print("Approximation: ", rounded_element)
                    print("Exact: ", end="")
                    display(element)
        else:
            print("Ellipse not in standard form")


class Hyperbola(Conic):
    def __init__(self, equation):
        super().__init__(equation)
        self.standard_form = False

    @property
    def centre(self):
        return self.get_centre()

    @property
    def orientation(self):
        """
        Return the orientation of the hyperbola: 'horizontal', 'vertical' or 'rotated'.
        """
        # Calculate the angle in radians
        theta = sympy.Rational(1, 2) * sympy.atan(sympy.Rational(self.B, (self.A - self.C)))

        # Use sympy's equivalence checking for comparison
        if self.B != 0:
            return 'Rotated', theta
        else:
            if self.A > self.C:
                return 'Horizontal', 0
            elif self.A < self.C:
                return 'Vertical', 0

    @property
    def semi_transverse_axis(self):
        """To be implemented: return the length of the semi-transverse axis"""
        pass

    @property
    def semi_conjugate_axis(self):
        """To be implemented: return the length of the semi-conjugate axis"""
        pass

    @property
    def eccentricity(self):
        """To be implemented: return the eccentricity of the hyperbola"""
        pass

    @property
    def vertices(self):
        """To be implemented: return the vertices of the hyperbola"""
        pass

    @property
    def foci(self):
        """To be implemented: return the foci of the hyperbola"""
        pass

    @property
    def directrices(self):
        """To be implemented: return the directrices of the hyperbola"""
        pass

    @property
    def asymptotes(self):
        """To be implemented: return the equations of the asymptotes of the hyperbola"""
        pass

    def get_info(self):
        print(f"{self.__repr__()}\nType: {self.type}\nCoefficients: {self.coeff}"
              f"\nGeneral Form: {self}\n")
        self.print_matrices()
        print(f"Centre: {self.centre}")
        if self.orientation[0] == 'rotated':
            angle_in_degrees = sympy.N(sympy.deg(self.orientation[1]), 5)
            print(f"Orientation: {self.orientation[0]}, {angle_in_degrees} degrees\nradians: {self.orientation[1]}\n")
        else:
            print(f"Orientation: {self.orientation[0]}\n")
        print(f"Equation of transverse axis")
        display(self.transverse_axis_line())
        self.plot()

    def get_centre(self):
        """
        The center of the hyperbola is the point where the derivatives w.r.t x and y are zero.
        It's a stationary point of the hyperbola
        :return: tuple(x_coordinate, y_coordinate)
        """
        # Define the symbols
        x, y = sympy.symbols('x y')

        # Define the system of equations

        # eq1: partial derivative w.r.t x
        eq1 = sympy.Eq(2 * self.A * x + self.B * y + self.D, 0)
        # eq2: partial derivative w.r.t y
        eq2 = sympy.Eq(self.B * x + 2 * self.C * y + self.E, 0)

        # Solve the system
        solution = sympy.solve((eq1, eq2), (x, y))

        center_x = solution[x]
        center_y = solution[y]

        return center_x, center_y

    def transverse_axis_line(self):
        """
        Return the equation of the line passing through the transverse axis of the hyperbola.
        """
        x, y = sympy.symbols('x y')

        # Calculate the slope from the rotation angle
        if self.orientation[0] == 'Rotated':
            m = sympy.tan(self.orientation[1])
            # Calculate the y-intercept from the y-coordinate of the center
            c = self.centre[1] - m * self.centre[0]
            eqn = sympy.nsimplify(sympy.Eq(y, m * x + c), 0.0001)
            simplified_eq = sympy.Eq(eqn.lhs, sympy.trigsimp(eqn.rhs))
            return simplified_eq

        elif self.orientation[0] == 'Horizontal':
            return sympy.Eq(y, self.centre[1])
        else:  # Vertical
            return sympy.Eq(x, self.centre[0])

    def _get_rotation_angle(self):
        """Abstract method of rotate"""
        pass

    def _get_translation_point(self):
        """Abstract method of translate"""
        h, k = self.centre[0], self.centre[1]
        original = h, k
        return h, k, original

    def _standard_form_flag(self):
        """Abstract method of rotate and translate"""
        # Update standard_form flag
        if self.centre == (0, 0) and self.orientation[0] == 'Horizontal':
            self.standard_form = True
            # Eradicate the constant term in coefficient matrix (Normalise the equation)
            self.coeff_matrix = abs(1/self.F) * self.coeff_matrix
            # Update coefficients
            self._update_from_matrix()

    def record_state(self):
        """To be implemented: record the current state of the hyperbola"""
        pass

    def print_history(self):
        """To be implemented: print the history of the transformations applied to the hyperbola"""
        pass

    def plot(self, x_range=None, y_range=None):
        """
        Method to plot an instance of a hyperbola along with its centre and transverse axis line.
        """
        plt.figure(figsize=(10, 8))

        # Default ranges
        if x_range is None:
            x_range = (-10, 10)
        if y_range is None:
            y_range = (-10, 10)

        x = np.linspace(*x_range, 400)
        y = np.linspace(*y_range, 400)
        x, y = np.meshgrid(x, y)

        # Plot the hyperbola
        plt.contour(x, y, (self.A * x ** 2 + self.B * x * y + self.C * y ** 2
                           + self.D * x + self.E * y + self.F), [1], colors='r')

        # Plot the centre
        plt.plot(*self.centre, 'go')

        # Plot the transverse axis line
        if self.orientation[0] == 'Rotated':
            transverse_axis_line = self.transverse_axis_line()
            x_line = np.linspace(*x_range, 100)
            y_line = sympy.lambdify(sympy.symbols('x'), transverse_axis_line.rhs, 'numpy')(x_line)
            plt.plot(x_line, y_line, 'b--')
        elif self.orientation[0] == 'Horizontal':
            plt.axhline(self.centre[1], color='b', linestyle='--')
        else:  # Vertical
            plt.axvline(self.centre[0], color='b', linestyle='--')

        plt.gca().set_aspect('equal', adjustable='box')
        plt.axhline(0, color='gray', linewidth=0.5)
        plt.axvline(0, color='gray', linewidth=0.5)
        plt.title(f"${self.__str__()}$")
        plt.show()


