# v11 manuscript placement and A4 print review

The current manuscript is Letter (215.9 x 279.4 mm), not A4. In the disposable
RevTeX preview, textwidth inserts the 6.75-in native figures at
179.25 mm, or 104.55% of native size.
The main glyphs are at least 2.56 mm;
math scripts are approximately 1.77 mm.
PNG effective resolution at that placement is 573.9 dpi.

| Figure | Page | Width (mm) | Height (mm) |
| --- | ---: | ---: | ---: |
| figure1_relaxation_times_comparison | 6 | 179.25 | 188.54 |
| figure2_transport_coefficients_comparison | 7 | 179.25 | 156.67 |

| A4 printing assumption | Scale | Figure width (mm) | Main glyph (mm) | Script (mm) |
| --- | ---: | ---: | ---: | ---: |
| A4_fit_page | 97.27% | 174.35 | 2.49 | 1.72 |
| A4_fit_5mm_printable_margin | 92.64% | 166.05 | 2.38 | 1.64 |

Fit-page is a geometric page-box calculation. A real printer may impose a
different printable area; the 5-mm case is an explicit example, not a claim
about the user's printer. Actual-size printing keeps the Letter PDF glyph
sizes, subject to the printer's printable-area limits.

The temporary source updates only paths, captions and figure float options.
The original [t] setting produced float-placement warnings for the taller
figure/caption; [!t] in the preview allows the available page area and keeps
10 pages without those warnings. The formal paper files are not changed.

Visual placement is usable at the tested width. This is not formal submission
compliance: math scripts remain below 2 mm, textwidth slightly exceeds the
repository 7-in profile cap, and an enlarged PNG falls below 600 effective
dpi. The review PDFs remove only the raster-resolution limitation, not the
font-size or width gates. Do not silently promote manuscript eligibility.

The two A4 print-preview PDFs are disposable, two-page excerpts of the vector
manuscript preview. See placement_report.json for absolute paths and hashes.
