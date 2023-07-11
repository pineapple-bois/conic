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
        x = np.linspace(-2, 2, 400)
        y = np.linspace(-2, 2, 400)
        x, y = np.meshgrid(x, y)

        plt.contour(x, y, (self.A * x ** 2 + self.B * x * y + self.C * y ** 2
                           + self.D * x + self.E * y + self.F), [1], colors='r')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.axhline(0, color='gray', linewidth=0.5)
        plt.axvline(0, color='gray', linewidth=0.5)
        plt.title(f"${self.__str__()}$")
        plt.show()


class Parabola(Conic):
    def __init__(self, equation):
        super().__init__(equation)

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
        print(f"\nAngle of rotation: {self.get_rotation_angle()}")
        print(f"Orientation: {self.get_orientation()}")
        print(f"Axis of symmetry: {self.axis}")
        self.plot()

    def get_orientation(self):
        # Determine the orientation based on the values of A, B, and C
        if self.A != 0 and self.C == 0:
            if self.A > 0:
                return "vertical", "positive"
            else:
                return "vertical", "negative"
        elif self.A == 0 and self.C != 0:
            if self.C > 0:
                return "horizontal", "positive"
            else:
                return "horizontal", "negative"
        elif self.A != 0 and self.C != 0:
            # For a rotated parabola, return 'R' for now.
            # You will need to compute the angle of rotation here.
            theta = 0.5 * sympy.atan2(self.B, (self.A - self.C))
            return "rotated", theta

    def get_rotation_angle(self):
        # TODO: METHOD - find out which quadrant based on signs of coefficients
        if self.B ** 2 - 4 * self.A * self.C != 0:
            raise ValueError('The conic section is not a parabola')

        # Avoid division by zero
        if self.B != 0:
            theta = sympy.Rational(0.5) * sympy.acot((sympy.Rational(self.A - self.C)) / sympy.Rational(self.B))
            # convert to numerical result
            theta = sympy.N(theta)
            degrees = sympy.N(theta * (180/sympy.pi), 3)
            return theta, degrees
        else:
            print("Not Rotated")

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
                h = sympy.Rational(-d / (2 * a)).limit_denominator(100000)
                k = sympy.Rational(a * h ** 2 + d * h + f).limit_denominator(100000)

            else:  # Parabola opens to the right or left
                c = self.C / -self.D
                e = self.E / -self.D
                f = self.F / -self.D
                k = sympy.Rational(-e / (2 * c)).limit_denominator(100000)
                h = sympy.Rational(c * k ** 2 + e * k + f).limit_denominator(100000)

        else:  # handle rotated case
            # The logic for computing vertex for rotated parabolas will go here.
            # This depends on how you plan to compute it. If you're using the axis property
            # and solve method as in the commented function, you can add those lines here.
            x, y = sympy.symbols('x y')
            gen_eqn = self.save_as_sympy(return_expr=True)
            axis_y = sympy.solve(self.axis, y)[0]  # replace 'axis' with the appropriate variable/method
            axis_eqn = sympy.Eq(y, axis_y)

            solution = sympy.solve((gen_eqn, axis_eqn), (x, y))
            if solution:
                h, k = solution[0]
            else:
                h, k = None, None
                print("No solution found.")

        return h, k

    def plot(self):
        h, k = self.vertex
        margin = 5  # adjust this to your needs

        if self.orientation[0] == 'vertical':  # Parabola opens up or down
            x = np.linspace(float(h) - margin, float(h) + margin, 400)
            y = (self.A * x ** 2 + self.D * x + self.F) / -self.E
            plt.plot(x, y, color='r')
            plt.axvline(x=self.axis, color='b', linestyle='dotted')

        elif self.orientation[0] == 'horizontal':  # Parabola opens to the right or left
            y = np.linspace(float(k) - margin, float(k) + margin, 400)
            x = (self.C * y ** 2 + self.E * y + self.F) / -self.D
            plt.plot(x, y, color='r')
            plt.axhline(y=self.axis, color='b', linestyle='dotted')

        elif self.orientation[0] == 'rotated':
            # Plot the parabola
            x = np.linspace(-20, 20, 400)
            y = np.linspace(-20, 20, 400)
            x, y = np.meshgrid(x, y)
            plt.contour(x, y, (self.A * x ** 2 + self.B * x * y + self.C * y ** 2
                            + self.D * x + self.E * y + self.F), [1], colors='r')

            # Plot the axis of symmetry
            x_sym, y_sym = sympy.symbols('x y')
            a = self.axis.coeff(x_sym)
            b = self.axis.coeff(y_sym)
            m = -a / b  # slope
            c = self.axis.subs({x_sym: 0, y_sym: 0}) / b  # intercept
            x_values = np.linspace(-20, 20, 400)
            y_values = m * x_values - c
            plt.plot(x_values, y_values, color='b', linestyle='dotted')

        plt.gca().set_aspect('equal', adjustable='box')
        plt.axhline(0, color='gray', linewidth=0.5)
        plt.axvline(0, color='gray', linewidth=0.5)
        plt.plot(h, k, 'bo')  # plot the vertex as a blue dot
        plt.annotate(f'Vertex', (h, k), textcoords="offset points", xytext=(-10, -10), ha='center')  # Annotate vertex
        plt.title(f"General Form:\n${self}$")
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



