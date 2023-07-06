import re


def decompose_polynomial(formula: str):
    """
    The input formula should be a string representation of the polynomial.

    Polynomial exponents are natural numbers.
    If exponents are rationals or negative integers, the formula will not parse

    :return Polynomial dictionary;
    {(exponent_x, exponent_y, exponent_z): coefficient}
    """

    if not isinstance(formula, str):
        raise TypeError("Input must be a string representation of a polynomial.")

    formula = formula.replace(' ', '').replace('^', '**')
    # Extract all variables in formula
    variables = set(re.findall('[a-z]', formula))

    # Check that the variables are only x, y, and z
    if variables.difference({'x', 'y', 'z'}):
        raise ValueError("For multivariate polynomials, variables can only be 'x', 'y', or 'z'")

    terms = re.findall(
        r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:\*[a-z](?:\*\*\d+)?)*(?:\*[a-z](?:\*\*\d+)?)*|[+-]?'
        r'(?:[a-z](?:\*\*\d+)?)*(?:\*[a-z](?:\*\*\d+)?)*|[+-]?[a-z]', formula)

    terms = [term for term in terms if term]  # remove any empty strings
    terms = [term.replace('**', '^') for term in terms]  # analysis is simpler with only one iteration of '*'
    # Constructing the dictionary
    poly_dict = {}
    for term in terms:
        if term.startswith('+'):
            term = term[1:]

        if '*' in term:
            parts = term.split('*')

            # Extract and handle the coefficient if it exists
            if parts[0].replace('.', '').isdigit() or parts[0].replace('-', '').replace('.', '').isdigit():
                coefficient = float(parts.pop(0))
            else:
                coefficient = -1.0 if parts[0].startswith('-') else 1.0
                if parts[0].startswith('-'):
                    parts[0] = parts[0][1:]

            # Initialize list of exponents
            exponents = [0, 0, 0]

            for part in parts:
                # Each part is either of the form 'x', 'x^2', etc.
                if '^' in part:
                    variable, exponent = part.split('^')
                    exponent = int(exponent)
                else:
                    variable = part
                    exponent = 1

                if variable == 'x':
                    exponents[0] = exponent
                elif variable == 'y':
                    exponents[1] = exponent
                elif variable == 'z':
                    exponents[2] = exponent

            # The key in the dictionary is a tuple of exponents
            tuple_key = tuple(exponents)
            poly_dict[tuple_key] = coefficient

        elif '^' in term:  # Handle terms like '+x^2', '+y^2', etc.
            parts = term.split('^')
            variable = parts[0]
            if variable.startswith('-'):
                variable = variable[1:]
                coefficient = -1.0
            else:
                coefficient = 1.0
            exponent = int(parts[1])

            exponents = [0, 0, 0]

            if variable == 'x':
                exponents[0] = exponent
            elif variable == 'y':
                exponents[1] = exponent
            elif variable == 'z':
                exponents[2] = exponent

            tuple_key = tuple(exponents)
            poly_dict[tuple_key] = coefficient

        elif term.isalpha() or (len(term) > 1 and term[1:].isalpha()):
            # statement for terms like '+x', '-y' etc ... no explicit coefficient
            variable = term[-1]
            coefficient = 1.0 if term[0] != '-' else -1.0
            exponents = [0, 0, 0]

            if variable == 'x':
                exponents[0] = 1
            elif variable == 'y':
                exponents[1] = 1
            elif variable == 'z':
                exponents[2] = 1

            tuple_key = tuple(exponents)
            if tuple_key not in poly_dict:
                poly_dict[tuple_key] = coefficient

        elif term.lstrip('+-').replace('.', '').isdigit():
            # statement for constant terms like '+2', '-1.5', '3' etc.
            coefficient = float(term)
            poly_dict[(0, 0, 0)] = coefficient

    return poly_dict

