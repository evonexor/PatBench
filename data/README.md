# data/

This directory is empty on purpose: **no data is stored in the repository.**

- The data workspace (patent archives and the output of every stage) and every release directory live outside the repository. Point `PATBENCH_WORKSPACE` and `--out` at them. The layout is described in [docs/REPRODUCE.md](../docs/REPRODUCE.md).
- `.gitignore` excludes `*.parquet`, `*.zip` and everything in this directory except this file.
- To run the pipeline on a small sample, put it here (git ignores it) and point `PATBENCH_WORKSPACE` at it.

What a release directory contains and the columns of every table: [docs/DATA_SCHEMA.md](../docs/DATA_SCHEMA.md).
Splits and task sets: [docs/SPLITS.md](../docs/SPLITS.md).
