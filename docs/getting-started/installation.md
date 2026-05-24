# Installation

## Requirements

- Python >= 3.10, < 3.13

## Install from Source

```bash
git clone https://github.com/5h4ng/OptiMHC.git
cd OptiMHC
pip install -e .
```

Or using **uv** (recommended for development):

```bash
git clone https://github.com/5h4ng/OptiMHC.git
cd OptiMHC
uv sync --locked
```

!!! note "PyPI"
    PyPI distribution is not yet available. For now, install from source as shown above.

## NetMHCpan / NetMHCIIpan Setup

Most OptiMHC dependencies are installed automatically via pip. The only exceptions are **NetMHCpan** and **NetMHCIIpan**, which are standalone executables that must be downloaded and installed manually.

### 1. Download

- [NetMHCpan 4.1](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/) (MHC Class I)
- [NetMHCIIpan 4.3](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/) (MHC Class II)

Both require a license from DTU Health Tech (free for academic use).

### 2. Install and Add to PATH

After downloading and extracting, follow the installation instructions provided by DTU. Then make sure the executables are accessible from your shell:

=== "Linux / macOS"

    ```bash
    # Add to your ~/.bashrc or ~/.zshrc:
    export PATH="/path/to/netMHCpan-4.1:$PATH"
    export PATH="/path/to/netMHCIIpan-4.3:$PATH"
    ```

=== "Verify"

    ```bash
    which netMHCpan       # Should print the path to the executable
    which netMHCIIpan     # Should print the path to the executable
    ```

OptiMHC calls these tools via the [mhctools](https://github.com/openvax/mhctools) library, which invokes `netMHCpan` and `netMHCIIpan` as command-line programs. If they are not on your `PATH`, the prediction step will fail with a "command not found" error.

You can also set the executable path per feature generator instead of modifying `PATH`:

```yaml
featureGenerator:
  - name: NetMHCpan
    params:
      executablePath: /path/to/netMHCpan
  - name: NetMHCIIpan
    params:
      executablePath: /path/to/netMHCIIpan
```

### 3. Test the Installation

```bash
netMHCpan -v          # Should print the version number
netMHCIIpan -v        # Should print the version number
```

!!! tip
    If you do not need MHC binding predictions, you can skip this step entirely and use OptiMHC with other features (Basic, SpectralSimilarity, DeepLC, PWM, OverlappingPeptide, MHCflurry) — all of which are installed via pip.

## Verify OptiMHC Installation

```bash
optimhc --help
```

You should see the available commands: `pipeline` and `experiment`.
