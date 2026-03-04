.PHONY: install lint test format docs build

install:
	poetry install

lint:
	poetry run ruff check src tests
	poetry run mypy src

format:
	poetry run black src tests
	poetry run ruff check src tests --fix

test:
	poetry run pytest --cov=mmalignments --cov-report=term-missing

docs:
	poetry run sphinx-build docs docs/_build/html

build:
	poetry build
