# Image prompts for the documentation

Generate each image, save it in this folder under the **File name** given, and say
when they are in. Nothing references them yet, the vignettes still carry the
ASCII diagrams, and each one is swapped for its image only once the file exists.

## House style, applies to every prompt below

Every prompt already carries this, but if you regenerate one, keep it:

- **Flat vector illustration, no photorealism, no 3D, no gradients, no drop
  shadows.** These sit beside body text and have to read at 700 px wide.
- **Palette:** background transparent or white. Primary accent `#EC7414`
  (the logo's orange). Secondary `#0072F0` blue, `#00A651` green, `#8E44E8`
  purple. Greys `#333333` for text and axes, `#CCCCCC` for gridlines. **Never
  red-and-green as the only distinction between two things**, around 8% of male
  readers cannot separate them.
- **Type:** clean sans-serif, large enough to read at 700 px. Label directly on
  the drawing rather than in a legend wherever possible.
- **Aspect:** 16:9 unless the prompt says otherwise.
- **No invented numbers.** Where a prompt asks for axis labels, use exactly the
  text given. An illustration that puts a plausible-looking figure on an axis is
  a figure a reader will quote.

A note on what these are for: they replace ASCII diagrams that already work. The
image has to be *clearer* than the ASCII, not merely prettier, if a generated
one is ambiguous where the ASCII was not, the ASCII stays.

---

## 1. `fluidics-single-file.png`

**Replaces:** the prose in *Flow cytometry for dummies* §1, which currently has no
diagram at all.

> Flat vector diagram, 16:9, white background. A narrow vertical fluid stream
> flowing downward through the centre. Inside the stream, six round cells spaced
> in single file, queueing one behind another. A horizontal laser beam in orange
> `#EC7414` crosses the stream at one point, and exactly one cell sits in the beam
> at that moment, drawn glowing. To the right, an arrow leads to a simple table
> with three rows and four columns, the row aligned with the lit cell highlighted
> in orange, showing that one crossing becomes one row of numbers. Label the
> stream "cells in single file", the beam "laser", the table "one row per event".
> Clean sans-serif labels in dark grey `#333333`. No 3D, no gradients.

---

## 2. `two-humps-and-the-cut.png`

**Replaces:** the ASCII histogram in §3, "The cut".

> Flat vector line chart, 16:9, white background. A single smooth curve showing a
> bimodal distribution: a tall left hump and a shorter right hump, with a clear
> dip between them. The area under the left hump filled pale grey `#CCCCCC`, under
> the right hump filled orange `#EC7414` at about 40% opacity. A vertical dashed
> dark grey line dropped exactly at the lowest point of the dip between the humps,
> labelled "the cut". The left hump labelled "negative, no marker", the right
> hump labelled "positive, marker present". X axis labelled "brightness", Y axis
> labelled "number of cells". No numbers on either axis. Clean sans-serif.

---

## 3. `valley-vs-fixed-cut.png`

**Replaces:** nothing yet. New, for §3, this is the idea the whole package rests
on and no diagram states it.

> Flat vector diagram, 16:9, white background, two panels side by side separated
> by a thin vertical grey rule. Each panel shows two bimodal curves stacked
> vertically, one labelled "sample A" and one labelled "sample B", where sample B's
> curves are shifted noticeably to the right of sample A's. Left panel titled "one
> fixed cut for every sample": a single vertical dashed line at the same x position
> in both curves, correctly splitting sample A but landing in the middle of sample
> B's right hump, with that error circled in `#0072F0` blue and labelled "lands
> inside the positive population". Right panel titled "each sample gets its own
> cut": two vertical dashed orange `#EC7414` lines, each sitting in its own curve's
> dip, both correct. Clean sans-serif, no numbers on axes.

---

## 4. `fsc-ssc-populations.png`

**Replaces:** the ASCII scatter in §6, "FSC and SSC".

> Flat vector scatter plot, 16:9, white background. X axis labelled "FSC, size",
> Y axis labelled "SSC, internal graininess", both with an arrow head and no
> numbers. Four labelled clusters of small dots: "debris" bottom-left, small and
> scattered, in light grey; "lymphocytes" low and left, tight cluster, blue
> `#0072F0`; "monocytes" middle, medium cluster, purple `#8E44E8`; "granulocytes"
> upper right, large cluster, orange `#EC7414`. Each cluster labelled directly
> beside it, not in a legend. A short annotation with a leader line pointing at the
> granulocyte cluster reading "grainiest, which is where the name comes from".
> Clean sans-serif.

---

## 5. `pulse-a-h-w.png`

**Replaces:** the pulse description in §7, "The `-A`, `-H` and `-W` suffixes".

> Flat vector diagram, 16:9, white background. A single smooth bell-shaped pulse
> curve drawn in dark grey, with time on the X axis labelled "time as the cell
> crosses the beam" and signal on the Y axis. Three measurements annotated on the
> same curve in three different colours: the area under the curve shaded orange
> `#EC7414` and labelled "-A, area"; a vertical line from the baseline to the peak
> in blue `#0072F0` labelled "-H, height"; a horizontal double-headed arrow across
> the base of the pulse in green `#00A651` labelled "-W, width". Clean sans-serif,
> no numbers on axes.

---

## 6. `doublet-vs-singlet.png`

**Replaces:** nothing yet. New, for §7, the reason `-W` exists is easier to see
than to read.

> Flat vector diagram, 16:9, white background, two panels side by side. Left panel
> labelled "one cell": a single round cell crossing a laser beam, and beneath it a
> narrow bell-shaped pulse, with the pulse width marked by a short green `#00A651`
> double-headed arrow labelled "narrow". Right panel labelled "two cells stuck
> together": two round cells touching, crossing the same beam, and beneath them a
> wider pulse with a slight double peak, its width marked by a longer green arrow
> labelled "wide, this is how a doublet is caught". Same axis treatment in both
> panels. Clean sans-serif, no numbers.

---

## 7. `gating-hierarchy.png`

**Replaces:** nothing yet. New, for §9, "Declare the big populations first".

> Flat vector funnel or nested-boxes diagram, 4:3, white background. Four stages
> narrowing downward, each a horizontal bar narrower than the one above it, in a
> single orange `#EC7414` at decreasing opacity from dark at the top to pale at the
> bottom. Top to bottom the bars are labelled "all events", "single cells",
> "live cells", "CD45 positive, white blood cells", and beneath the last a
> branching row of three small boxes labelled "T cells", "B cells", "monocytes".
> To the right of the funnel, a vertical annotation with a downward arrow reading
> "an error at any level is inherited by everything below it". Clean sans-serif.

---

## 8. `compositional-constraint.png`

**Replaces:** nothing yet. New, for the *Statistics* article, section 3
"Compositionality". This is the single most misread idea in the output and it has
no picture.

> Flat vector diagram, 16:9, white background. Two horizontal stacked bars of
> identical total length, one above the other, each divided into four coloured
> segments labelled inside: "granulocytes" orange `#EC7414`, "T cells" blue
> `#0072F0`, "B cells" green `#00A651`, "NK cells" purple `#8E44E8`. The top bar
> labelled "healthy" with roughly equal segments. The bottom bar labelled "patient"
> where the granulocyte segment is much larger and all three others are
> proportionally squeezed. Downward arrows from each shrunken segment in the top
> bar to its counterpart in the bottom bar, each marked with a small downward
> triangle. A caption box to the right reading "only granulocytes changed. The
> other three fell because the bar is a fixed length." Clean sans-serif, both bars
> exactly the same width.

---

## 9. `donors-not-cells.png`

**Replaces:** nothing yet. New, for *Statistics* §1, "Pseudoreplication".

> Flat vector diagram, 16:9, white background, two panels side by side. Left panel
> titled "counting cells": six person icons in dark grey, each with a large stack
> of many small dots beneath them, and a total beneath reading "n = 1,600,000",
> with a red-free warning marker, use a purple `#8E44E8` cross rather than red,
> and the label "wrong: n is set by how long the machine ran". Right panel titled
> "counting donors": the same six person icons, each with a single orange
> `#EC7414` dot beneath, total reading "n = 6", with a green `#00A651` tick and the
> label "right: n is the number of people". Clean sans-serif.

---

## 10. `interval-not-point.png`

**Replaces:** nothing yet. New, for *Statistics* §4b, on clinical effect intervals.

> Flat vector diagram, 16:9, white background. A horizontal axis running from -1 on
> the left to +1 on the right, with 0 marked by a vertical dashed grey line
> labelled "no effect". Three rows stacked vertically, each with a filled circle at
> its estimate and a horizontal bar through it for the interval. Row 1 labelled
> "n = 9": circle at about +0.6 with a very wide bar running from about -0.15 to
> +1.0, drawn in grey, annotated "a lead, not a result". Row 2 labelled "n = 40":
> circle at the same +0.6 with a much narrower bar from about +0.35 to +0.8, drawn
> in orange `#EC7414`, annotated "the same effect, now pinned down". Row 3 labelled
> "n = 9" with a circle at +0.05 and an equally wide grey bar, annotated
> "indistinguishable from the row above, on this cohort". Clean sans-serif.

---
---

# Part 2. Concepts explained elsewhere on the site

Part 1 covers *Flow cytometry for dummies* and *Statistics*. The prompts below
cover every other place the documentation explains a cytometry concept in prose or
ASCII and would read better as a picture. Same house style throughout.

---

## 11. `transform-arcsinh-vs-logicle.png`

**For:** *Gating specification* section 4, Transformation. The choice between the
two is explained in words and the difference is entirely visual.

> Flat vector line chart, 16:9, white background. Two curves on shared axes
> mapping raw detector value on the X axis to transformed value on the Y axis.
> The X axis crosses zero, running from a negative region on the left through zero
> into a large positive region, labelled "raw channel value", with a pale grey
> shaded band over the negative region labelled "negative values, real after
> compensation". Curve 1 in orange EC7414 labelled "arcsinh": smooth, symmetric
> through zero, near-linear close to zero and compressing at the extremes. Curve 2
> in blue 0072F0 labelled "logicle": similar at the positive end but with a
> visibly wider near-linear region around zero. An annotation with a leader line
> into the grey band reading "both keep negatives; a plain log cannot". Clean
> sans-serif, no numbers on axes.

---

## 12. `gate-hierarchy-inheritance.png`

**For:** *Gating specification* section 1, Hierarchy.

> Flat vector tree diagram, 4:3, white background. A vertical parent-child tree.
> Root node "all events" at the top in dark grey. Below it "scatter gate", below
> that "singlets", below that "live", below that "CD45 positive". From the CD45
> node three children branch out: "T cells", "B cells", "myeloid". One node in the
> chain, "live", is outlined in orange EC7414 with a thick border and a small
> warning glyph. Every node beneath it is drawn faded to 40 percent opacity. An
> annotation to the right with a leader line to the outlined node reads "get this
> cut wrong and every number below it is wrong, with no error raised". Clean
> sans-serif.

---

## 13. `threshold-drift.png`

**For:** *Diagnostics* section 4, Threshold drift. The failure mode is subtle and
a picture makes it obvious.

> Flat vector dot plot, 16:9, white background. X axis has two categories,
> "healthy" and "patient". Y axis labelled "where the cut landed", no numbers.
> Under healthy, six orange EC7414 dots clustered low. Under patient, six blue
> 0072F0 dots clustered clearly higher, with no overlap between the clusters. A
> horizontal bracket spanning the gap annotated "the cuts themselves differ by
> group". Below the plot a caption box reading "then part of any difference in the
> result is a difference in the definition, not the biology". Clean sans-serif.

---

## 14. `gate-uncertainty-resample.png`

**For:** *Diagnostics* section 5, Gate uncertainty.

> Flat vector diagram, 16:9, white background, three parts left to right. Left: one
> bimodal density curve with a single solid orange EC7414 vertical line in its dip,
> labelled "the cut this run reported". Centre: a right-pointing arrow labelled
> "resample the cells and re-derive, many times". Right: the same curve with about
> twenty faint vertical grey lines scattered across a narrow band around the dip,
> the band shaded and marked with a horizontal double-headed arrow labelled "how
> far the cut can move". A caption beneath reading "the second number reported
> beside every frequency". Clean sans-serif, no numbers.

---

## 15. `lod-loq.png`

**For:** *Diagnostics* section 5a, Detection limits.

> Flat vector diagram, 16:9, white background. A horizontal axis labelled
> "population frequency", increasing to the right, no numbers. Three zones divided
> by two vertical dashed lines. Leftmost zone shaded dark grey, labelled "below LOD
>, cannot be told apart from nothing". Middle zone shaded pale orange, labelled
> "detected but not quantifiable, present, but the number is not reliable". Right
> zone white, labelled "quantified". The two dividing lines labelled "LOD" and
> "LOQ". Beneath, an annotation reading "a rare population reported as 0.02 percent
> may be any of these three". Clean sans-serif.

---

## 16. `batch-vs-group-confounding.png`

**For:** *Diagnostics* section 10, Batch structure. The one diagnostic that can
refuse to run, and the reason is geometric.

> Flat vector diagram, 16:9, white background, two panels side by side. Each panel
> is a 2 by 2 grid of cells, rows labelled "batch 1" and "batch 2", columns
> labelled "healthy" and "patient". Left panel titled "separable": all four cells
> contain three small sample icons, evenly filled, annotated beneath "correction is
> possible, each group appears in each batch". Right panel titled "confounded":
> only the top-left and bottom-right cells contain icons, the other two empty,
> annotated beneath "correction is refused, removing the batch and removing the
> finding are the same operation". Orange EC7414 icons for healthy, blue 0072F0 for
> patient. Clean sans-serif.

---

## 17. `cluster-gate-agreement.png`

**For:** *Diagnostics* section 6, Cluster concordance.

> Flat vector diagram, 16:9, white background. A scatter of cells arranged as four
> rounded blobs, with an overlay of dashed outlines representing declared gates.
> Blob 1 fully enclosed by one dashed outline, labelled "declared and found,
> agree". Blob 2 split across two dashed outlines, labelled "one declared label
> spanning two clusters". Blob 3 containing two separate dashed outlines, labelled
> "two labels inside one cluster". Blob 4 with no dashed outline at all, filled
> orange EC7414, labelled "found, never declared, what explore mode is for".
> Clean sans-serif, no axis numbers.

---

## 18. `gate-transferability.png`

**For:** *Diagnostics* section 7, Gate transferability.

> Flat vector diagram, 16:9, white background, two panels. Left panel titled
> "held-out cells": six donor icons all orange EC7414 inside a dashed box labelled
> "fitted here", with an arrow to a score reading "F1 0.94" and the annotation "the
> same donors, flattering". Right panel titled "held-out donor": five orange donor
> icons inside the dashed "fitted here" box and one blue 0072F0 donor icon outside
> it, arrow to a score reading "F1 0.71", annotation "a donor the gate never saw,
> the honest number". Clean sans-serif.

---

## 19. `pipeline-ten-stages.png`

**For:** *How it works, the ten stages*. One map of the whole pipeline.

> Flat vector flow diagram, 16:9, white background. Ten numbered boxes in two rows
> of five, snaking so the flow reads row one left to right then row two left to
> right, with arrows between consecutive boxes. In order: 1 Read, 2 Transform,
> 3 Gate, 4 Score, 5 Uncertainty, 6 Embed, 7 Test, 8 Diagnose, 9 Explain,
> 10 Report. Boxes 1 to 5 orange EC7414, box 6 purple 8E44E8, box 7 blue 0072F0,
> boxes 8 to 10 green 00A651. A detached box below, connected by a dashed line,
> labelled "10a Explore, optional, changes nothing above it". A small annotation
> under stage 3 reading "per sample, from its own density". Clean sans-serif.

---

## 20. `two-input-files.png`

**For:** *Inputs, the sample sheet and config*. The split between the two files is
what new users get wrong.

> Flat vector diagram, 16:9, white background. Two document icons side by side.
> Left document titled "samples.csv" with an orange EC7414 header, showing a small
> table with the column headers "file", "sample_id", "patient_id", "cohort",
> "batch" and three rows beneath. Caption under it: "one row per FCS file,
> anything that VARIES per sample". Right document titled "analysis.yaml" with a
> blue 0072F0 header, showing indented key-value lines suggesting "populations:",
> "gating:", "ratios:". Caption under it: "one decision for the whole study". A
> vertical dashed divider between them with a note reading "an analysis choice in
> the CSV would be repeated on every row and invite the rows to disagree". Clean
> sans-serif.

---

## 21. `attainable-p-floor.png`

**For:** *Every output file*, on group_differences.png. The most counter-intuitive
statement in the documentation, and it needs a picture.

> Flat vector scatter plot styled as a volcano plot, 16:9, white background. X axis
> labelled "effect (Cliff's delta)" running from -1 to +1 with 0 marked. Y axis
> labelled "-log10(p)". About twelve grey dots scattered across the lower half. A
> red dashed horizontal line low down labelled "p = 0.05". A blue 0072F0 dot-dash
> horizontal line drawn ABOVE the red one, labelled "the smallest p this design can
> reach". The band between the two lines shaded pale grey with text across it
> reading "no population can land here, however cleanly the groups separate". A
> caption beneath reading "an empty top half is a fact about the design, not about
> the biology". Clean sans-serif.

---

## 22. `explore-vs-declared.png`

**For:** *Explore mode, unsupervised discovery*.

> Flat vector diagram, 16:9, white background, two side-by-side panels sharing an
> identical background scatter of cells. Left panel titled "declared": four regions
> outlined and filled orange EC7414 at low opacity, each labelled with a population
> name, and a large remaining region left plain grey labelled "82 percent of cells
>, no name". Right panel titled "explore": the same scatter divided instead into
> about twelve irregular coloured cluster regions covering everything including the
> previously grey area, three of them ringed in purple 8E44E8 and labelled "no
> declared population covers these". An arrow between the panels labelled "ignores
> the specification and the parent gate". Clean sans-serif.

---

## 23. `per-sample-vs-per-patient.png`

**For:** *Statistics* section 4a, the unit of analysis for clinical variables. This
is the defect that produced two false findings on a real cohort, so it earns a
diagram.

> Flat vector diagram, 16:9, white background, two panels. Both show the same three
> person icons, each with three test-tube icons beneath it, nine tubes in total.
> Left panel titled "SOFA, measured at each draw": each of the nine tubes carries
> a different small number badge in orange EC7414, annotated "the value moves
> within a patient, so the sample is the unit, n = 9". Right panel titled
> "survival, a property of the person": the three tubes under each person all
> carry the SAME badge as each other in blue 0072F0, the three people differing
> from one another, annotated "counting nine would be counting three people three
> times each, n = 3". Clean sans-serif.

---

## 24. `spillover-spreading.png`

**For:** *Gating specification* section 5, Compensation, and the spreading report.

> Flat vector diagram, 16:9, white background. Two panels stacked vertically
> sharing an X axis labelled "detector B brightness". Top panel titled "cells
> negative for marker A": a single narrow bell curve centred near zero, its width
> marked by a short double-headed arrow labelled "narrow", drawn in grey. Bottom
> panel titled "cells bright for marker A": a bell curve at the same centre
> position but visibly wider, width marked by a longer arrow labelled "wider, this
> is spreading", drawn in orange EC7414. An annotation to the right reading
> "compensation puts the centre back. It cannot put the width back." Clean
> sans-serif, no numbers.

---

## 25. `event-vs-cell.png`

**For:** *Flow cytometry for dummies* section 1, on why the word is "event".

> Flat vector diagram, 16:9, white background. Three panels in a row, each showing
> one thing crossing a laser beam and each producing one row in a small table to
> the right. Panel 1: a single clean round cell, labelled "one cell, what you
> want". Panel 2: two cells stuck together, labelled "a doublet, counted as one
> event". Panel 3: an irregular fragment, labelled "debris from a dead cell". All
> three rows in the table are drawn identically and marked with a single brace
> labelled "the file cannot tell them apart, the gates have to". Cells in orange
> EC7414, debris in grey. Clean sans-serif.
