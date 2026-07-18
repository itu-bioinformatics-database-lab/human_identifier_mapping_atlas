#!/usr/bin/env bash
# Push human_identifier_mapping_atlas to GitHub.
# Usage:
#   bash push_to_github.sh YOUR_GITHUB_USERNAME [REPO_NAME]
# Example:
#   bash push_to_github.sh itu-bioinformatics-database-lab human_identifier_mapping_atlas

set -euo pipefail

USER_OR_ORG="${1:?Usage: $0 <github-username-or-org> [repo-name]}"
REPO_NAME="${2:-human_identifier_mapping_atlas}"
PKG_DIR="$(cd "$(dirname "$0")" && pwd)"

cd "$PKG_DIR"

if [[ ! -d .git ]]; then
  git init
  git branch -M main
fi

git add -A
git status

echo
echo "Large raw/mapping files are gitignored (see DATA_SOURCES.md)."
read -r -p "Commit message [Initial commit: human identifier mapping atlas]: " MSG || true
MSG="${MSG:-Initial commit: human identifier mapping atlas}"

if git diff --cached --quiet; then
  echo "Nothing new to commit."
else
  git commit -m "$MSG"
fi

REMOTE_URL="https://github.com/${USER_OR_ORG}/${REPO_NAME}.git"

if git remote get-url origin >/dev/null 2>&1; then
  echo "Remote origin: $(git remote get-url origin)"
  git push -u origin main
  echo "Done: https://github.com/${USER_OR_ORG}/${REPO_NAME}"
  exit 0
fi

if command -v gh >/dev/null 2>&1; then
  gh repo create "${USER_OR_ORG}/${REPO_NAME}" --public --source=. --remote=origin --push \
    --description "Human gene/protein ID atlas: Ensembl ↔ gene symbol ↔ UniProt mappings (JSON), with regenerable BioMart/HGNC/UniProt downloads."
  echo "Done: https://github.com/${USER_OR_ORG}/${REPO_NAME}"
  exit 0
fi

git remote add origin "$REMOTE_URL"
echo
echo "Remote set to: $REMOTE_URL"
echo "Create an empty repo on GitHub named '${REPO_NAME}' (if needed), then run:"
echo "  cd \"$PKG_DIR\""
echo "  git push -u origin main"
