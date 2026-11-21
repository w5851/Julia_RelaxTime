# Frozen v15 transport backgrounds

These point tables are the published inputs used by the reviewed manuscript
`paper/My_Paper/pnjl_aniso_transport`. They contain 942 phase-referenced rows
and 925 fixed-temperature rows, including the original branch labels and
the closest sampled points at first-order transitions.

The original numerical calculation used commit
`2e5c411c716d5cbc524c3fc6cfb335fe05b39444`. Source intake:
`writing_package/attachments/data_sources/v2_onshellkernel/publication_clean_v15_pdf/intake_manifest.json`.
`input_manifest.json` binds both original CSV byte hashes and their scan groups.

The Ward/Møller rerun reuses the masses, Polyakov loops, densities, entropy,
temperature, chemical potential, anisotropy and branch identities from each
row. Original relaxation times and transport values are retained only for
terminal-kernel replay and comparison. They are not results of the new
prescription. New outputs use a separate result case.
