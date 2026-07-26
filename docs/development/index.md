# Development

This page covers how to set up a development environment for contributing to OptiMHC.

## Clone the Repository

```bash
git clone https://github.com/5h4ng/OptiMHC.git
cd OptiMHC
```

## Install Dependencies

OptiMHC uses [uv](https://docs.astral.sh/uv/) for dependency management (recommended):

```bash
uv sync --locked --group dev
```

Alternatively, using pip:

```bash
pip install -e .
pip install pytest ruff pre-commit
```

## Running Tests

```bash
uv run pytest                              # All tests
uv run pytest tests/test_psm_container_contract.py  # Single file
uv run pytest tests/ -k "test_config"      # Filter by name
```

## Linting and Formatting

OptiMHC uses [Ruff](https://docs.astral.sh/ruff/) for linting and formatting:

```bash
uv run ruff check .          # Lint
uv run ruff format .         # Format
uv run ruff format --check . # Check without modifying
```

Configuration: `line-length = 99`, rules `["E", "F", "I"]` (pycodestyle errors,
pyflakes, isort). `E501` (line too long) is ignored. Ruff excludes `docs/`.

## Pre-commit Hooks

```bash
uv run pre-commit install              # Install hooks
uv run pre-commit run --all-files      # Run all hooks manually
```

## Building Documentation

```bash
mkdocs serve   # Local preview at http://127.0.0.1:8000
mkdocs build   # Static build to site/
```
