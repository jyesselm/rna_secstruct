.PHONY: help install install-dev test test-cov lint format type-check clean build publish publish-test

help:
	@echo "Available commands:"
	@echo "  make install       - Install package"
	@echo "  make install-dev   - Install package with dev dependencies"
	@echo "  make test          - Run tests"
	@echo "  make test-cov      - Run tests with coverage"
	@echo "  make lint          - Run linter (ruff)"
	@echo "  make format        - Format code (black)"
	@echo "  make format-check  - Check code formatting"
	@echo "  make type-check    - Run type checker (mypy)"
	@echo "  make check-all     - Run all checks (lint, format, type, test)"
	@echo "  make clean         - Remove build artifacts"
	@echo "  make build         - Build package"
	@echo "  make publish       - Build and publish to PyPI"
	@echo "  make publish-test  - Build and publish to Test PyPI"

install:
	pip install -e .

install-dev:
	pip install -e ".[dev,all]"

test:
	pytest

test-cov:
	pytest --cov=rna_secstruct --cov-report=html --cov-report=term --cov-report=xml

test-fast:
	pytest -m "not slow"

lint:
	ruff check rna_secstruct/ test/

lint-fix:
	ruff check --fix rna_secstruct/ test/

format:
	black rna_secstruct/ test/

format-check:
	black --check rna_secstruct/ test/

type-check:
	mypy rna_secstruct/

check-all: format-check lint type-check test
	@echo "All checks passed!"

clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info
	rm -rf .pytest_cache/
	rm -rf .mypy_cache/
	rm -rf .ruff_cache/
	rm -rf htmlcov/
	rm -rf .coverage
	rm -rf coverage.xml
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

build: clean
	python -m build

publish: check-all build
	@echo "Publishing to PyPI..."
	python -m twine upload dist/*

publish-test: check-all build
	@echo "Publishing to Test PyPI..."
	python -m twine upload --repository testpypi dist/*

