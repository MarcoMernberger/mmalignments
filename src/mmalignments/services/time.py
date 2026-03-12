"""A module to take care of time, PID and formatting."""

from datetime import datetime

# The format string for timestamps.
TIMEFORMAT = "%Y-%m-%d-%H-%M-%S"


def timestamp() -> datetime:
    """
    Get the current timestamp.

    Returns
    -------
    str
        _description_
    """
    return datetime.now()


def timestamp_to_str(timestamp: datetime) -> str:
    """
    Convert a datetime object to a string format used in log filenames.

    Parameters
    ----------
    timestamp : datetime
        The timestamp to convert.

    Returns
    -------
    str
        The formatted timestamp string.
    """
    return timestamp.strftime(TIMEFORMAT)


def str_to_timestamp(cur_ts_str: str) -> datetime:
    """
    Convert a timestamp string from log filenames back to a datetime object.

    Parameters
    ----------
    cur_ts_str : str
        The timestamp string to convert.

    Returns
    -------
    datetime
        The corresponding datetime object.
    """
    ret = datetime.strptime(cur_ts_str, TIMEFORMAT)
    return ret
