import numpy as np
import matplotlib.pyplot as plt
import sympy
from fractions import Fraction
from IPython.display import display
from matplotlib.lines import Line2D
from poly_dictionary import decompose_polynomial


class Conic:
    def __init__(self, equation: str):
        self.original_input = equation
        self.equation_str, self.equation = self._to_general_form_conic(equation)
        self.coefficients = decompose_polynomial(self.equation_str)
        self.coefficients = self._remove_z_coordinate(self.coefficients)
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
        equation = str(equation).replace('**', '^').replace('*', '').replace('sqrt', '√')
        return equation

    def create_coeff_matrix(self):
        """
        Create the matrix of coefficients,
        """
        matrix = np.array([[self.A, self.B / 2, self.D / 2],
                           [self.B / 2, self.C, self.E / 2],
                           [self.D / 2, self.E / 2, self.F]])
        return matrix

    def create_quad_matrix(self):
        """
        Create 2x2 sub-matrix from coefficient matrix
        """
        matrix = np.array([[self.A, self.B / 2],
                           [self.B / 2, self.C]])
        return matrix

    def print_matrices(self):
        print(f"Matrix of coefficients:\n{self.coeff_matrix}\n\nQuadratic Matrix:\n{self.quad_matrix}")

    def classify(self):
        """
        Classify the conic section based on its coefficients.
        Using the determinant of the 2x2 sub-matrix (delta) and the determinant of the 3x3 matrix (Delta).
        """
        # Calculate the determinant of the 2x2 quadratic-matrix
        delta = np.linalg.det(self.quad_matrix)

        # TODO: Methods remain, useful for classification however, unused
        # Calculate trace of the 2x2 sub-matrix
        tau = np.trace(self.quad_matrix)

        # Calculate the determinant of the 3x3 matrix. - DEGENERATE CONICS?
        Delta = np.linalg.det(self.coeff_matrix)

        if delta == 0:
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

    def _to_general_form_conic(self, equation: str):
        """
        Takes the string, works out the LCM, collects like terms
        then reduces the resultant as far as possible.
        :return: General form equation of a conic as string and sympy expression
        """
        x, y = sympy.symbols('x y')
        equation_str = sympy.sympify(equation, locals={"sqrt": sympy.sqrt})

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
        formula = sympy.expand(formula)     # General form requires ZERO arbitrary factorisation
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

    def draw(self):
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

    # TODO: NEED TO DEAL WITH SPECIAL ANGLES IN [30, 45, 60] DEGREES
    # TODO: NEED TO UPDATE THE INFO AND EQUATION WITHIN SUBCLASS
    # TODO: NEED TO WORK ON THE TRANSLATION TO STANDARD POSITION & WRITE IN FORM Y^2=4AX

class Parabola(Conic):
    def __init__(self, equation):
        super().__init__(equation)
        self.history = []

    @property
    def orientation(self):
        return self.get_orientation()

    @property
    def axis(self):
        return self.compute_axis()

    @property
    def vertex(self):
        return self.compute_vertex()

    def get_info(self):
        print(f"{self.__repr__()}\nType: {self.type}\nCoefficients: {self.coeff}"
              f"\nGeneral Form: {self}\n")
        self.print_matrices()
        print(f"\nOrientation: {self.get_orientation()}")
        print(f"Axis of symmetry: {self.axis}")
        self.plot()

    def get_orientation(self):
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

    def rotate_parabola(self, rational=False):
        orientation, rotation_angle = self.get_orientation()

        if orientation == "vertical":
            if rotation_angle == "positive":
                rotation_angle = sympy.pi * 3 / 2  # 270 degrees in radians
            elif rotation_angle == "negative":
                rotation_angle = sympy.pi / 2  # 90 degrees in radians
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
        self.update_from_matrix()
        self.record_state()

        # After rotating, if the parabola isn't oriented as 'horizontal, positive',
        # recursively rotate again
        orientation, rotation_angle = self.get_orientation()
        if orientation == "vertical" or (orientation == "horizontal" and rotation_angle == "negative"):
            self.rotate_parabola(rational)

    def update_coefficients(self):
        self.A = self.coeff_matrix[0, 0]
        self.B = 2 * self.coeff_matrix[0, 1]  # Remember we stored B/2 in the matrix
        self.C = self.coeff_matrix[1, 1]
        self.D = 2 * self.coeff_matrix[0, 2]  # Remember we stored D/2 in the matrix
        self.E = 2 * self.coeff_matrix[1, 2]  # Remember we stored E/2 in the matrix
        self.F = self.coeff_matrix[2, 2]

    def update_from_matrix(self):
        # Update coefficients from the matrix
        self.update_coefficients()

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
        # Store the current state in history
        self.history.append({
            'Coefficients': self.coefficients,
            'Matrix': self.coeff_matrix,
            'Equation': self.equation,
            'Orientation': self.orientation
        })

    def print_history(self):
        print(f"Original input:\n{self.original_input}\n")
        for item in self.history:
            if isinstance(item, str):
                print(f"----------\n{item}\n")
            elif isinstance(item, dict):
                for key, value in item.items():
                    print(f"{key} :\n{value}\n")

    def compute_axis(self):
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
            alpha = np.sqrt(self.A)
            if self.B >= 0:
                beta = np.sqrt(self.C)
            else:
                beta = -np.sqrt(self.C)
            axis = alpha * gen_x + beta * gen_y

        return axis

    def compute_vertex(self):
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

    def plot(self, x_range=None, y_range=None):
        # Increase figure size
        plt.figure(figsize=(10, 8))  # Adjust the values as per your desired size

        if self.orientation[0] == 'rotated':
            if not x_range:
                x_range = [-10, 10]
            if not y_range:
                y_range = [-10, 10]
        else:
            h, k = self.vertex
            if not x_range:  # use default range if not provided
                x_range = [float(h) - 5, float(h) + 5]
            if not y_range:
                y_range = [float(k) - 5, float(k) + 5]

        if self.orientation[0] == 'vertical':  # Parabola opens up or down
            x = np.linspace(x_range[0], x_range[1], 400)
            y = (self.A * x ** 2 + self.D * x + self.F) / -self.E
            plt.plot(x, y, color='r')
            plt.axvline(x=self.axis, color='b', linestyle='dotted')
            plt.plot(h, k, 'bo')  # plot the vertex as a blue dot
            plt.annotate(f'Vertex', (h, k), textcoords="offset points",
                         xytext=(-10, -10), ha='center')  # Annotate vertex
            custom_lines = [Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10),
                            Line2D([0], [0], color='blue', linestyle='dotted')]
            axis_label = f"x={self.axis}"
            plt.legend(custom_lines, [f'Vertex: ({h}, {k})', axis_label], loc='best')


        elif self.orientation[0] == 'horizontal':  # Parabola opens to the right or left
            y = np.linspace(y_range[0], y_range[1], 400)
            x = (self.C * y ** 2 + self.E * y + self.F) / -self.D
            plt.plot(x, y, color='r')
            plt.axhline(y=self.axis, color='b', linestyle='dotted')
            plt.plot(h, k, 'bo')  # plot the vertex as a blue dot
            plt.annotate(f'Vertex', (h, k), textcoords="offset points",
                         xytext=(-10, -10), ha='center')  # Annotate vertex
            custom_lines = [Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10),
                            Line2D([0], [0], color='blue', linestyle='dotted')]
            axis_label = f"y={self.axis}"
            plt.legend(custom_lines, [f'Vertex: ({h}, {k})', axis_label], loc='best')

        elif self.orientation[0] == 'rotated':
            # This will require both x and y ranges
            x = np.linspace(x_range[0], x_range[1], 400)
            y = np.linspace(y_range[0], y_range[1], 400)
            x, y = np.meshgrid(x, y)
            plt.contour(x, y, (self.A * x ** 2 + self.B * x * y + self.C * y ** 2
                               + self.D * x + self.E * y + self.F), [1], colors='r')

            # Plot the axis of symmetry
            x_sym, y_sym = sympy.symbols('x y')
            a = self.axis.coeff(x_sym)
            b = self.axis.coeff(y_sym)
            m = -a / b  # slope
            c = self.axis.subs({x_sym: 0, y_sym: 0}) / b  # intercept
            x_values = np.linspace(x_range[0], x_range[1], 400)
            y_values = m * x_values - c
            plt.plot(x_values, y_values, color='b', linestyle='dotted')

        plt.gca().set_aspect('equal', adjustable='box')
        plt.xlim(x_range)
        plt.ylim(y_range)
        plt.axhline(0, color='gray', linewidth=0.5)
        plt.axvline(0, color='gray', linewidth=0.5)
        plt.title(f"${str(self)}$")
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



