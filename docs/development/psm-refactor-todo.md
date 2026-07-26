# PSM refactor follow-ups

- [ ] Move the legacy automatic pepXML p-value/E-value log heuristic out of `read_pepxml` into an explicit, independently tested `SearchScoreTransform`. Keep the first bug-fix refactor behavior-compatible and do not expand it with this extraction.
- [ ] Upgrade Mokapot after it publishes a release containing the `config.seed` forwarding now present on its `main` branch (commit `5bad097`, PR #143). OptiMHC is pinned to Mokapot 0.10.0, whose CLI sets only NumPy's legacy global seed and does not pass `rng` to the model or `brew`.
