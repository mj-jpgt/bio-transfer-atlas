env:
	mamba env create -f environment.yml

download:
	python scripts/download_1000g.py
	python scripts/download_pgs.py
	python scripts/download_reactome.py

smoke:
	snakemake --snakefile workflow/Snakefile \
	  --configfile configs/experiments/mvp_chr22.yaml \
	  --cores 4

full:
	snakemake --snakefile workflow/Snakefile \
	  --configfile configs/experiments/full_autosomes.yaml \
	  --cores 16

figures:
	snakemake results/figures/atlas_main.pdf --cores 4

test:
	pytest -q

manifest:
	python scripts/update_manifest.py

clean_interim:
	rm -rf data/interim/

.PHONY: env download smoke full figures test manifest clean_interim
