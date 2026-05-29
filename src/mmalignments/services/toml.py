

def indent_toml(text: str, indent: str = "    ") -> str:
    """
    Indents the lines of a TOML document.

    Parameters
    ----------
    text : str
        TOML document as a string.
    indent : str, optional
        String used for indentation, by default "    "

    Returns
    -------
    str
        Indented TOML document.
    """
    lines = text.splitlines()

    result = []
    inside_table = False

    for line in lines:
        stripped = line.strip()

        # table header
        if stripped.startswith("[") and stripped.endswith("]"):
            inside_table = True
            result.append(line)
            continue

        # leerzeile beendet block
        if stripped == "":
            inside_table = False
            result.append(line)
            continue

        # kommentare innerhalb table ebenfalls einrücken
        if inside_table:
            result.append(indent + line)
        else:
            result.append(line)

    return "\n".join(result)


def indent_regular_tables(
    text: str,
    indent: str = "    ",
) -> str:
    """
    Indents the lines of a TOML document that are within regular tables 
    (not array-of-tables).

    Parameters
    ----------
    text : str
        TOML document as a string.
    indent : str, optional
        String used for indentation, by default "    "

    Returns
    -------
    str
        Indented TOML document.
    """
    lines = text.splitlines()

    result = []
    inside_regular_table = False

    for line in lines:
        stripped = line.strip()

        # array-of-table => NICHT einrücken
        if stripped.startswith("[[") and stripped.endswith("]]"):
            inside_regular_table = False
            result.append(line)
            continue

        # normale table
        if stripped.startswith("[") and stripped.endswith("]"):
            inside_regular_table = True
            result.append(line)
            continue

        # leerzeile beendet block
        if stripped == "":
            inside_regular_table = False
            result.append(line)
            continue

        if inside_regular_table:
            result.append(indent + line)
        else:
            result.append(line)

    return "\n".join(result)