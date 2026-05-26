try:
    import string

    maketrans = string.maketrans
except (ImportError, NameError, AttributeError):
    maketrans = bytes.maketrans


rev_comp_table = maketrans(
    b"ACBDGHKMNSRUTWVYacbdghkmnsrutwvy", b"TGVHCDMKNSYAAWBRTGVHCDMKNSYAAWBR"
)


def reverse_complement(sequence: str) -> str:
    """
    reverse_complement returns the reverse complement of given sequence.

    Parameters
    ----------
    sequence : str
        Input sequence.

    Returns
    -------
    str
        Reverse complement of input sequence.
    """
    return sequence[::-1].translate(rev_comp_table)
