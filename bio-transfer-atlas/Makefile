env:
	mamba env create -f environment.yml

download:
	python scripts/download_1000g.py
	python scripts/download_pgs.py
	python scripts/download_reactome.py

test:
	pytest -q

manifest:
	python scripts/update_manifest.py

clean_interim:
	rm -rf data/interim/

.PHONY: env download test manifest clean_interim
