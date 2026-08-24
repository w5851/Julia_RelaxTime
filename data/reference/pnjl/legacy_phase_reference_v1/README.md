# Legacy PNJL phase-reference snapshot v1

This directory is the byte-preserving retirement snapshot for the historical
Issue #130 phase-reference inputs. The files are no longer canonical
repository-root inputs; they remain available for explicit rollback and for
per-key fallback when the strict candidate has no certified row.

The files keep their historical schemas. `RETIREMENT_MANIFEST.json` records
the pre-retirement SHA-256 and byte count. No solver or numerical regeneration
was used for this move.
