import tarfile
from pathlib import Path

dest = Path("/lambda/nfs/geeg/fairness")
with tarfile.open("/tmp/bta_scripts_update.tgz") as t:
    for m in t.getmembers():
        t.extract(m, path=dest, set_attrs=False)
print("scripts_synced")
