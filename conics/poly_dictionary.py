import re


def decompose_polynomial(formula: str):
    """
    The input formula should be a string representation of the polynomial.

    Conic exponents are natural numbers.
    If exponents are rationals or negative integers, the formula will not parse

    :return Polynomial dictionary;
    {(exponent_x, exponent_y): coefficient}
    """
    if not isinstance(formula, str):
        raise TypeError("Input must be a string representation of a conic.")

    formula = formula.replace(' ', '').replace('^', '**')
    variables = set(re.findall('[a-z]', formula))

    if variables.difference({'x', 'y'}):
        raise ValueError("Variables can only be 'x' or 'y'")

    terms = re.findall(
        r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:\*[a-z](?:\*\*\d+)?)*|[+-]?'
        r'(?:[a-z](?:\*\*\d+)?)*|[+-]?[a-z]', formula)

    terms = [term for term in terms if term]  # remove any empty strings
    terms = [term.replace('**', '^') for term in terms]

    poly_dict = {}
    for term in terms:
        if term.startswith('+'):
            term = term[1:]

        if '*' in term:
            parts = term.split('*')

            if parts[0].replace('.', '').isdigit() or parts[0].replace('-', '').replace('.', '').isdigit():
                coefficient = int(parts.pop(0))
            else:
                coefficient = -1 if parts[0].startswith('-') else 1
                if parts[0].startswith('-'):
                    parts[0] = parts[0][1:]

            exponents = [0, 0]

            for part in parts:
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

            tuple_key = tuple(exponents)
            poly_dict[tuple_key] = coefficient

        elif '^' in term:
            parts = term.split('^')
            variable = parts[0]
            if variable.startswith('-'):
                variable = variable[1:]
                coefficient = -1
            else:
                coefficient = 1
            exponent = int(parts[1])

            exponents = [0, 0]

            if variable == 'x':
                exponents[0] = exponent
            elif variable == 'y':
                exponents[1] = exponent

            tuple_key = tuple(exponents)
            poly_dict[tuple_key] = coefficient

        elif term.isalpha() or (len(term) > 1 and term[1:].isalpha()):
            variable = term[-1]
            coefficient = 1 if term[0] != '-' else -1
            exponents = [0, 0]

            if variable == 'x':
                exponents[0] = 1
            elif variable == 'y':
                exponents[1] = 1

            tuple_key = tuple(exponents)
            if tuple_key not in poly_dict:
                poly_dict[tuple_key] = coefficient

        elif term.lstrip('+-').replace('.', '').isdigit():
            coefficient = int(term)
            poly_dict[(0, 0)] = coefficient

    return poly_dict

