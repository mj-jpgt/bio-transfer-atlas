# Create GitHub repo and push (run from fairness/ after authenticating).
# Usage:
#   $env:Path = "C:\Users\mobar\AppData\Local\Temp\gh-cli\bin;$env:Path"
#   gh auth login
#   powershell -File scripts/push_to_github.ps1

$ErrorActionPreference = "Stop"
$gh = (Get-Command gh -ErrorAction SilentlyContinue)?.Source
if (-not $gh) {
  $cand = "C:\Users\mobar\AppData\Local\Temp\gh-cli\bin\gh.exe"
  if (Test-Path $cand) { $gh = $cand } else { throw "gh not found; install GitHub CLI or add it to PATH" }
}

& $gh auth status
$login = & $gh api user --jq .login
$repo = "bio-transfer-atlas"
$name = "$login/$repo"

# Create if missing
$exists = & $gh repo view $name 2>$null
if (-not $exists) {
  & $gh repo create $repo --public --source=. --remote=origin --description "Biological Transferability Atlas: open-data PRS/GWAS portability risk + honest interventions"
} else {
  git remote remove origin 2>$null
  git remote add origin "https://github.com/$name.git"
}

git push -u origin HEAD
Write-Host "Pushed: https://github.com/$name"
