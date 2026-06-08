#!/usr/bin/env bash
set -euo pipefail

OLD_YEAR=2023
TARGET_YEAR=2026
HOLDER="Oak Ridge National Laboratory"

OLD_NOTICE="Copyright (C) ${OLD_YEAR} ${HOLDER}"
NEW_NOTICE="Copyright (C) ${OLD_YEAR}-${TARGET_YEAR} ${HOLDER}"

APPLY=0

if [[ "${1:-}" == "--apply" ]]; then
    APPLY=1
elif [[ "${1:-}" == "" || "${1:-}" == "--dry-run" ]]; then
    APPLY=0
else
    echo "Usage: $0 [--dry-run|--apply]" >&2
    exit 2
fi

git rev-parse --is-inside-work-tree >/dev/null

BASE_REF="${BASE_REF:-}"

if [[ -z "$BASE_REF" ]]; then
    BASE_REF="$(
        git rev-list -n 1 --before='2024-01-01 00:00:00' HEAD || true
    )"
fi

if [[ -z "$BASE_REF" ]]; then
    echo "Could not determine a baseline commit before 2024-01-01." >&2
    echo "Set BASE_REF manually, for example:" >&2
    echo "  BASE_REF=<tag-or-commit> $0 --dry-run" >&2
    exit 2
fi

first_added_year() {
    local file="$1"

    git log \
        --follow \
        --find-renames \
        --diff-filter=A \
        --format='%ad' \
        --date=format:%Y \
        -- "$file" |
    tail -n 1
}

has_nontrivial_changes_since_base() {
    local file="$1"

    ! git diff \
        --quiet \
        -w \
        --ignore-blank-lines \
        -I 'Copyright \(C\) [0-9,-]* Oak Ridge National Laboratory' \
        -I 'Version:[[:space:]]*[0-9]+' \
        -I '^[[:space:]]*-+[[:space:]]*$' \
        "$BASE_REF"..HEAD \
        -- "$file"
}

printf 'Using baseline: %s\n\n' "$BASE_REF"
printf '%-8s %-18s %s\n' "YEAR" "ACTION" "FILE"
printf '%-8s %-18s %s\n' "----" "------" "----"

while IFS= read -r -d '' file; do
    if ! has_nontrivial_changes_since_base "$file"; then
        printf '%-8s %-18s %s\n' "-" "skip-no-change" "$file"
        continue
    fi

    first_year="$(first_added_year "$file")"

    if [[ -z "$first_year" ]]; then
        printf '%-8s %-18s %s\n' "unknown" "skip-review" "$file"
        continue
    fi

    if (( first_year > OLD_YEAR )); then
        printf '%-8s %-18s %s -> %s\n' \
            "$first_year" "review-later-add" "$file" "$NEW_NOTICE"
        continue
    fi

    if (( APPLY )); then
        OLD_NOTICE="$OLD_NOTICE" NEW_NOTICE="$NEW_NOTICE" perl -0pi -e '
            my $old = quotemeta($ENV{"OLD_NOTICE"});
            my $new = $ENV{"NEW_NOTICE"};
            s/$old/$new/g;
        ' -- "$file"

        printf '%-8s %-18s %s -> %s\n' \
            "$first_year" "updated" "$file" "$NEW_NOTICE"
    else
        printf '%-8s %-18s %s -> %s\n' \
            "$first_year" "would-update" "$file" "$NEW_NOTICE"
    fi

done < <(git grep -Ilz --fixed-strings "$OLD_NOTICE" -- .)
