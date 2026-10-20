# BrownianMotion 0.4.0

## Breaking changes

* All exported functions now use snake_case names, so that no function name can
  be mistaken for an S3 method (e.g. `create.bm` previously read as a `create()`
  method for class `bm`). The old names have been removed:

  | Old name            | New name            |
  |---------------------|---------------------|
  | `create.bm()`       | `create_bm()`       |
  | `add.points()`      | `add_points()`      |
  | `add.segment()`     | `add_segment()`     |
  | `add.labels()`      | `add_labels()`      |
  | `delete.labels()`   | `delete_labels()`   |
  | `delete.skeleton()` | `delete_skeleton()` |
  | `concat.bm()`       | `concat_bm()`       |
  | `first.passage()`   | `first_passage()`   |

  `sim()`, `layers()`, `refine()` and `plot()` are unchanged.

## Bug fixes

* Adding localised layers before existing localised bridge intervals now
  preserves the existing intervals' transformed layer bounds.
* Simulating or refining within localised layers requested via
  `layers(type = "localised")` no longer fails with "missing value where
  TRUE/FALSE needed".
* Indexing a `BrownianMotion` object with exactly one layer now works both at
  and off skeleton points.
* Conditional simulation within intersection layers no longer fails due to a
  missing internal helper.
* Updated ggplot2 usage for compatibility with current ggplot2 releases, and
  added missing namespace imports.

## Internal

* Internal helper functions renamed to avoid clashing with S3 dispatch
  (e.g. `sim.condbessel_()` is now `sim_condbessel_()`, and `update.lyr_()`,
  which R treated as a `stats::update()` method, is now `update_lyr_()`).

# BrownianMotion 0.3.0

* New `add.labels()` and `delete.labels()` for managing user labels on existing
  path times; `delete.skeleton()` gains the ability to delete labels.
* Tutorial vignette completed.

# BrownianMotion 0.2.0

* Multi-dimensional Brownian motion with arbitrary covariance matrix, via
  `create.bm(dim = , cov = )`, with access to the transformed path through `$Z`.
* Unified layer construction through `layers()` (Bessel, intersection and
  localised layers).
* New `add.points()`, `add.segment()` and `concat.bm()`.
* `delete.skeleton()` for removing path observations and layers.
* Labelling of path times.
* `[` indexing of Brownian motion objects by time, index, time range and label.
* Multi-dimensional pairwise plotting (requires patchwork; `sf` suggested for
  union layer plots).

# BrownianMotion 0.1.0

* Initial version: unconditional simulation, Brownian bridges, first passage
  simulation and Bessel layers.
