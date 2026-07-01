.PHONY: lint fmt test test-py test-rs build clean pre-commit docs

lint: lint-py lint-rs  ## Run all linters

lint-py:
	cd python && ruff check .

lint-rs:
	cd rust && cargo clippy --workspace --all-targets -- -D warnings

fmt: fmt-py fmt-rs  ## Format all code

fmt-py:
	cd python && ruff format .

fmt-rs:
	cd rust && cargo fmt --all

test: test-py test-rs  ## Run all tests

test-py:
	cd python && pytest tests/ -v --tb=short

test-rs:
	cd rust && cargo test --workspace --all-targets

build: build-py build-rs  ## Build all packages

build-py:
	cd python && pip install -e .

build-rs:
	cd rust && cargo build --workspace

pre-commit:  ## Run pre-commit on all files
	pre-commit run --all-files

docs:  ## Build documentation locally
	mkdocs build

clean:  ## Remove build artifacts
	cd rust && cargo clean
	find python -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	rm -rf site/
