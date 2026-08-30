.PHONY: test lint validate dry-run

test:
	pytest

lint:
	ruff check src tests

validate:
	hostbias validate --config config/config.example.yaml

dry-run:
	snakemake --snakefile workflow/Snakefile --configfile config/config.example.yaml --dry-run all
