# Environment Configuration Notes

## FusionCatcher Environment Issues

The FusionCatcher conda environment requires special configuration to resolve dependency conflicts:

### Issue
FusionCatcher depends on legacy packages (openpyxl 2.5.0a2, Python 2.7) that have conflicts with modern conda repositories.

### Solution
When creating FusionCatcher environments, use flexible channel priority:

```bash
# For manual environment creation
mamba env create -f envs/fusioncatcher.yaml --channel-priority flexible

# For Snakemake (configure in ~/.mambarc or conda config)
mamba config --set channel_priority flexible
```

### Working Environment Configuration
The original `envs/fusioncatcher.yaml` works correctly:
```yaml
channels:
  - defaults
  - bioconda
  - conda-forge
dependencies:
  - python=2.7
  - fusioncatcher
```

### Performance Notes
- FusionCatcher database download: 4.4GB (completed successfully)
- Environment creation with flexible priority: ~93 packages, 714MB download
- All required dependencies resolve correctly with this configuration