import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import patches
import sympy

# PARABOLAS


def plot_parabola(conic, x_range=None, y_range=None):
    plt.figure(figsize=(10, 8))  # Adjust the values as per your desired size

    if conic.orientation[0] == 'Rotated':
        if not x_range:
            x_range = [-10, 10]
        if not y_range:
            y_range = [-10, 10]
    else:
        h, k = conic.vertex
        if not x_range:  # use default range if not provided
            x_range = [float(h) - 5, float(h) + 5]
        if not y_range:
            y_range = [float(k) - 5, float(k) + 5]

    if conic.orientation[0] == 'Vertical':  # Parabola opens up or down
        x = np.linspace(x_range[0], x_range[1], 400)
        y = (conic.A * x ** 2 + conic.D * x + conic.F) / -conic.E
        plt.plot(x, y, color='r')
        plt.axvline(x=conic.axis.rhs, color='b', linestyle='dotted')
        plt.plot(h, k, 'go')  # plot the vertex as a green dot
        custom_lines = [Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10),
                        Line2D([0], [0], color='blue', linestyle='dotted')]
        axis_label = f"${sympy.latex(conic.axis)}$"
        plt.legend(custom_lines, [f'Vertex: $({h}, {k})$', axis_label], loc='best')


    elif conic.orientation[0] == 'Horizontal':  # Parabola opens to the right or left
        y = np.linspace(y_range[0], y_range[1], 400)
        x = (conic.C * y ** 2 + conic.E * y + conic.F) / -conic.D
        plt.plot(x, y, color='r')
        plt.axhline(y=conic.axis.rhs, color='b', linestyle='dotted')
        plt.plot(h, k, 'go')  # plot the vertex as a green dot
        custom_lines = [Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10),
                        Line2D([0], [0], color='blue', linestyle='dotted')]
        axis_label = f"${sympy.latex(conic.axis)}$"
        plt.legend(custom_lines, [f'Vertex: $({h}, {k})$', axis_label], loc='best')

    elif conic.orientation[0] == 'Rotated':
        x = np.linspace(x_range[0], x_range[1], 400)
        y = np.linspace(y_range[0], y_range[1], 400)
        x, y = np.meshgrid(x, y)
        plt.contour(x, y, (conic.A * x ** 2 + conic.B * x * y + conic.C * y ** 2
                           + conic.D * x + conic.E * y + conic.F), [1], colors='r')

        # Plot the axis of symmetry
        x_sym, y_sym = sympy.symbols('x y')
        m = conic.axis.rhs.coeff(x_sym)
        b = conic.axis.rhs.subs(x_sym, 0)
        x_values = np.linspace(x_range[0], x_range[1], 400)
        y_values = m * x_values + b
        plt.plot(x_values, y_values, color='b', linestyle='dotted')

        custom_lines = [Line2D([0], [0], color='blue', linestyle='dotted')]
        axis_label = f"${sympy.latex(conic.axis)}$"
        plt.legend(custom_lines, [axis_label], loc='best')

    plt.gca().set_aspect('equal', adjustable='box')
    plt.xlim(x_range)
    plt.ylim(y_range)
    plt.axhline(0, color='gray', linewidth=0.5)
    plt.axvline(0, color='gray', linewidth=0.5)
    plt.title(f"General Form:\n${str(conic)}$")
    plt.show()


def parabola_standard(conic, x_range=None, y_range=None):
    plt.figure(figsize=(10, 8))  # Adjust the values as per your desired size
    if not conic.standard_form:
        raise ValueError("The parabola is not in standard form")

    if not x_range:
        x_range = [-5, 5]
    if not y_range:
        y_range = [-5, 5]

    y = np.linspace(y_range[0], y_range[1], 400)
    x = y**2 / (4 * conic.a)

    # Plot parabola
    plt.plot(x, y, label='Parabola')

    # Plot the radius of curvature
    theta = np.linspace(0, 2*np.pi, 1000)
    x_circle = 2*conic.a + 2*conic.a * np.cos(theta)
    y_circle = 2*conic.a * np.sin(theta)
    plt.plot(x_circle, y_circle, label='Radius of curvature')

    plt.scatter(conic.a, 0, color='r', label=f'Focus: ${conic.focus}$')  # plot the focus
    plt.plot([conic.a, conic.a], [-2*conic.a, 2*conic.a], color='g', label=f'Latus Rectum: ${conic.latus_rectum}$')  # plot the latus rectum
    plt.axvline(-conic.a, color='r', linestyle='--', label=f'Directrix: ${conic.directrix}$')  # plot the directrix
    plt.axhline(0, color='gray', linewidth=0.5)
    plt.axvline(0, color='gray', linewidth=0.5)
    plt.xlim(x_range)
    plt.ylim(y_range)
    plt.title(f'Standard Form:\n${conic}$')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.legend(loc='best')  # This will place the legend at the location that covers the least amount of data.
    plt.show()


# CIRCLES


def plot_circle(conic, x_range=None, y_range=None):
    plt.figure(figsize=(10, 8))  # Adjust the values as per your desired size
    if not x_range:
        x_range = [-10, 10]
    if not y_range:
        y_range = [-10, 10]

    theta = np.linspace(0, 2*np.pi, 1000)
    r = conic.radius

    x_centre, y_centre = conic.centre

    x = x_centre + r * np.cos(theta)
    y = y_centre + r * np.sin(theta)

    plt.plot(x, y, label='Circle')
    plt.scatter(x_centre, y_centre, color='r', label=f'Center: ({x_centre}, {y_centre})')  # plot the center
    plt.axhline(0, color='gray', linewidth=0.5)
    plt.axvline(0, color='gray', linewidth=0.5)
    plt.xlim(x_range)
    plt.ylim(y_range)
    plt.title(f'Standard Form:\n${conic}$')
    plt.legend(loc='best')  # This will place the legend at the location that covers the least amount of data.
    plt.gca().set_aspect('equal')  # Ensure aspect ratio is equal
    plt.show()


# ELLIPSES


def plot_ellipse(self, x_range=None, y_range=None):
    """
    Plot the ellipse with its center and semi-major axis.
    """
    # Increase figure size
    plt.figure(figsize=(10, 8))  # Adjust the values as per your desired size

    if not x_range:
        x_range = [-10, 10]
    if not y_range:
        y_range = [-10, 10]

    # Convert sympy expressions to lambdified expressions for plotting
    x, y = sympy.symbols('x y')
    lambda_axis = sympy.lambdify((x, y), self.semimajor_axis_line().lhs - self.semimajor_axis_line().rhs, 'numpy')

    x_vals, y_vals = np.mgrid[x_range[0]:x_range[1]:200j, y_range[0]:y_range[1]:200j]  # grid of points

    # Define the ellipse's equation from the coefficients
    z = (self.A * x_vals ** 2 + self.B * x_vals * y_vals + self.C * y_vals ** 2
         + self.D * x_vals + self.E * y_vals + self.F)

    # Plot the ellipse
    plt.contour(x_vals, y_vals, z, levels=[0], colors='r')

    # Plot the center
    plt.plot(self.centre[0], self.centre[1], 'go')

    # Plot the semi-major axis
    plt.contour(x_vals, y_vals, lambda_axis(x_vals, y_vals), levels=[0], colors='b',
                linewidths=0.5, linestyles='dashed')

    plt.gca().set_aspect('equal', adjustable='box')
    plt.axhline(0, color='gray', linewidth=0.5)
    plt.axvline(0, color='gray', linewidth=0.5)
    plt.title(f"${self.__str__()}$")
    plt.show()


def ellipse_standard(self, x_range=None, y_range=None, fig_size=(12, 8)):
    """
    Plot the ellipse with its foci, directrices, and vertices.
    """
    # Define the plot
    fig, ax = plt.subplots(figsize=fig_size)

    if not x_range:
        x_range = [-10, 10]
    if not y_range:
        y_range = [-10, 10]

    x_vals, y_vals = np.mgrid[x_range[0]:x_range[1]:200j, y_range[0]:y_range[1]:200j]  # grid of points

    # Define the ellipse's equation from the coefficients
    z = (self.A * x_vals ** 2 + self.C * y_vals ** 2 + self.F)

    # Plot the ellipse
    ax.contour(x_vals, y_vals, z, levels=[0], colors='red')

    # Plot the foci
    for focus in self.foci:
        ax.plot(focus[0], focus[1], 'go')

    # Plot the directrices
    for d in self.directrices:
        ax.axvline(d.rhs, color='blue', linewidth=0.5, linestyle='dashed')

    # Plot the vertices
    for v_set in self.vertices:
        for v in v_set:
            ax.plot(v[0], v[1], 'mo')

    ax.set_aspect('equal', adjustable='box')
    plt.axhline(0, color='gray', linewidth=0.5)
    plt.axvline(0, color='gray', linewidth=0.5)

    # Define the patches for the legend
    red_patch = patches.Patch(color='red', label='Ellipse')
    green_patch = patches.Patch(color='green', label='Foci')
    blue_patch = patches.Patch(color='blue', label=f'Directrices: $x= \\pm {sympy.latex(self.directrices[1].rhs)}$')
    magenta_patch = patches.Patch(color='magenta', label='Vertices')

    legend1 = plt.legend(handles=[red_patch, green_patch, magenta_patch,],
                         loc='upper center', bbox_to_anchor=(0.5, -0.10), ncol=3)
    plt.legend(handles=[blue_patch], loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2)
    plt.gca().add_artist(legend1)
    plt.subplots_adjust(left=0.1, bottom=0.2, right=0.9, top=0.9)
    plt.title(f"${self.__str__()}$")
    plt.tight_layout()  # To ensure that all elements fit within the figure boundaries
    plt.show()


# HYPERBOLAS

def plot_hyperbola(self, x_range=None, y_range=None, fig_size=(12, 8)):
    """
    Method to plot an instance of a hyperbola along with its centre and transverse axis line.
    """
    plt.figure(figsize=fig_size)

    # Default ranges
    if x_range is None:
        x_range = (-10, 10)
    if y_range is None:
        y_range = (-10, 10)

    # Convert sympy expressions to lambdified expressions for plotting
    x, y = sympy.symbols('x y')
    lambda_axis = sympy.lambdify((x, y), self.transverse_axis_line().lhs - self.transverse_axis_line().rhs, 'numpy')

    x_vals, y_vals = np.mgrid[x_range[0]:x_range[1]:200j, y_range[0]:y_range[1]:200j]  # grid of points

    # Define the hyperbola equation from the coefficients
    z = (self.A * x_vals ** 2 + self.B * x_vals * y_vals + self.C * y_vals ** 2
         + self.D * x_vals + self.E * y_vals + self.F)

    # Plot the hyperbola
    plt.contour(x_vals, y_vals, z, levels=[0], colors='r')

    # Plot the centre
    plt.plot(*self.centre, 'go')

    # Add legend
    custom_lines = [Line2D([0], [0], color='red', lw=2),
                    Line2D([0], [0], marker='o', color='w', markerfacecolor='green', markersize=10)]
    plt.legend(custom_lines, ['Hyperbola', f'Centre: $({self.centre[0]}, {self.centre[1]})$'], loc='best')

    # Plot the transverse axis line
    # We don't plot if in standard form as this line is the x-axis
    if not self.standard_form:

        if self.orientation[0] == 'Rotated':
            plt.contour(x_vals, y_vals, lambda_axis(x_vals, y_vals), levels=[0],
                        colors='b', linewidths=0.5, linestyles='dashed')

        elif self.orientation[0] == 'Horizontal':
            plt.axhline(self.centre[1], color='b', linestyle='--')
            # Update legend to include axis line
            custom_lines.append(Line2D([0], [0], color='blue', linestyle='--'))
            axis_label = f"${sympy.latex(self.transverse_axis_line())}$"
            plt.legend(custom_lines, ['Hyperbola', f'Centre: $({self.centre[0]},'
                                                   f' {self.centre[1]})$', axis_label], loc='best')

        else:  # Vertical
            plt.axvline(self.centre[0], color='b', linestyle='--')
            # Update legend to include axis line
            custom_lines.append(Line2D([0], [0], color='blue', linestyle='--'))
            axis_label = f"${sympy.latex(self.transverse_axis_line())}$"
            plt.legend(custom_lines, ['Hyperbola', f'Centre: $({self.centre[0]},'
                                                   f' {self.centre[1]})$', axis_label], loc='best')

    plt.gca().set_aspect('equal', adjustable='box')
    plt.axhline(0, color='gray', linewidth=0.5)
    plt.axvline(0, color='gray', linewidth=0.5)
    plt.title(f"${self.__str__()}$")
    plt.show()


def hyperbola_standard(self, x_range=None, y_range=None, fig_size=(16, 8)):
    """
    Plot the hyperbola with its foci, directrices, asymptotes, and vertices.
    """
    # Define the plot
    fig, ax = plt.subplots(figsize=fig_size)

    if not x_range:
        x_range = [-10, 10]
    if not y_range:
        y_range = [-10, 10]

    # Convert sympy expressions to lambdified expressions for plotting
    x, y = sympy.symbols('x y')

    x_vals, y_vals = np.mgrid[x_range[0]:x_range[1]:200j, y_range[0]:y_range[1]:200j]  # grid of points

    # Define the ellipse's equation from the coefficients
    z = (self.A * x_vals ** 2 + self.C * y_vals ** 2 + self.F)

    # Plot the hyperbola
    ax.contour(x_vals, y_vals, z, levels=[0], colors='red')

    # Plot the foci
    for focus in self.foci:
        ax.plot(focus[0], focus[1], 'go')

    # Plot the directrices
    for d in self.directrices:
        ax.axvline(d.rhs, color='blue', linewidth=0.5, linestyle='dashed')

    # Plot the asymptotes
    x_vals_asymptotes = np.linspace(x_range[0], x_range[1], 200)
    for a in self.asymptotes:
        y_vals_asymptotes = sympy.lambdify(x, a.rhs)(x_vals_asymptotes)
        ax.plot(x_vals_asymptotes, y_vals_asymptotes, color='purple', linewidth=0.5, linestyle='dashed')

    # Plot the vertices
    for vertex in self.vertices:
        ax.plot(vertex[0], vertex[1], 'co')

    ax.set_aspect('equal', adjustable='box')
    plt.axhline(0, color='gray', linewidth=0.5)
    plt.axvline(0, color='gray', linewidth=0.5)

    # Define the patches for the legend
    red_patch = patches.Patch(color='red', label='Hyperbola')
    green_patch = patches.Patch(color='green', label='Foci')
    blue_patch = patches.Patch(color='blue', label=f'Directrices: $x=\\pm {sympy.latex(self.directrices[1].rhs)}$')
    purple_patch = patches.Patch(color='purple', label=f'Asymptotes: $y=\\pm {sympy.latex(self.asymptotes[1].rhs)}$')
    cyan_patch = patches.Patch(color='cyan', label='Vertices')

    legend1 = plt.legend(handles=[red_patch, green_patch, cyan_patch,], loc='upper center', bbox_to_anchor=(0.5, -0.10), ncol=3)
    plt.legend(handles=[purple_patch, blue_patch], loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2)
    plt.gca().add_artist(legend1)
    plt.subplots_adjust(left=0.1, bottom=0.2, right=0.9, top=0.9)
    plt.title(f"${self.__str__()}$")
    plt.tight_layout()  # To ensure that all elements fit within the figure boundaries
    plt.show()
