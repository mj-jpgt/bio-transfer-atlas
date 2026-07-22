import tarfile
from pathlib import Path

root = Path("/lambda/nfs/geeg/fairness")
# score pgens -> interim
interim = root / "data/interim/1000g_grch38"
interim.mkdir(parents=True, exist_ok=True)
with tarfile.open("/tmp/bta_score_pgens.tgz") as t:
    for m in t.getmembers():
        t.extract(m, path=interim, set_attrs=False)
print("score_pgens_extracted", len(list(interim.glob("chr*.score.pgen"))))

# pgs + keeps at repo root
with tarfile.open("/tmp/bta_pgs_keeps.tgz") as t:
    for m in t.getmembers():
        t.extract(m, path=root, set_attrs=False)
print("pgs_keeps_extracted")
