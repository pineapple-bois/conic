# Conic Section Analyzer

----

This code provides a foundational structure for classifying and analyzing conic sections given an equation. Currently, the primary features include:

----

### Input Parsing and Processing: 
The code takes an equation as a string input, converts it to a Sympy expression, and reduces it to a polynomial form where the coefficients correspond to the general form of a conic section.

----

### Coefficient Extraction: 
The script extracts the coefficients of the equation and organizes them into a dictionary. The dictionary represents a simplified polynomial in the format coefficients = {(degree of x, degree of y): coefficient, ...}.

----

### Conic Classification: 
Based on the discriminant derived from these coefficients, the code classifies the given equation into one of four categories: Circle, Ellipse, Parabola, and Hyperbola.

----

### Object-oriented Design: 
The codebase utilizes an object-oriented design with a parent Conic class and subclasses for Circle, Ellipse, Parabola, and Hyperbola. This structure allows easy expansion for specific properties and methods applicable to each type of conic section.

----
### Visualization (in progress): 

The foundation for visualizing conic sections with Matplotlib has been laid, with the plot() method built into the base Conic class.
The main Conic class serves to extract key properties from the input equation and provide basic functions such as plotting the equation. Each subclass is meant to represent a specific type of conic section, and thus can be extended to incorporate properties and behaviors unique to that type.

----

## Future Work
The goal moving forward is to enhance the capability of this conic section analyzer by adding more features:

----

Transformation to Standard Form: Currently in progress, this feature will compute the appropriate linear transformations to convert the equation into its standard form. This process involves computing the transformation matrix and applying it to the given equation.

Subclass Specific Methods: Implementing methods such as plot() in each subclass (Parabola, Ellipse, Circle, and Hyperbola), to account for their specific characteristics when visualizing these conic sections.

Computational Properties: Additional computational properties unique to each type of conic section will be added in the subclasses. This could include details like the center, foci, and axes of an ellipse or the directrix and focus of a parabola.

-----

This conic section analyzer will serve as a robust tool for classifying, analyzing, and visualizing conic sections.
