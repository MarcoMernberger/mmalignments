import logging
from logging import Logger
from pathlib import Path

from .io import ensure
from .time import now_as_str


def initlog(console: bool = False, log_dir=None) -> Logger:
    timestamp = now_as_str()
    return setup_run_logger(timestamp, console=console, log_dir=log_dir)


def setup_run_logger(
    timestamp: str, console: bool = False, log_dir: Path | None = None
) -> Logger:
    log_dir = log_dir or Path(".logs")
    ensure(log_dir)
    run_log = log_dir / f"run_{timestamp}.log"
    logger = logging.getLogger("pipeline")
    logger.setLevel(logging.INFO)
    logger.propagate = False

    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s")

    if not any(
        isinstance(h, logging.FileHandler) and h.baseFilename == str(run_log)
        for h in logger.handlers
    ):
        fh = logging.FileHandler(run_log)
        fh.setLevel(logging.INFO)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    if console and not any(
        isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler)
        for h in logger.handlers
    ):
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(formatter)
        logger.addHandler(ch)

    return logger
