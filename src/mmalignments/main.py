import typer

from .config import settings
from .core.example import add

app = typer.Typer(help="Clean Python structure example")


@app.command()
def info() -> None:
    """Show basic environment info."""
    typer.echo(f"Environment: {settings.environment}")
    typer.echo(f"DB URL: {settings.db_url!r}")


@app.command()
def sum(a: int, b: int) -> None:
    """Add two integers."""
    result = add(a, b)
    typer.echo(f"{a} + {b} = {result}")


def main() -> None:
    app()


if __name__ == "__main__":
    main()
