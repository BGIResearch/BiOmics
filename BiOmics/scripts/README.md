# BiOmics-BRICK

BiOmics-BRICK is a lightweight command-line toolkit for the BRICK knowledge graph, providing querying, ranking/enrichment, graph conversion, visualization, and embedding. All functionality is available through the unified `brick` entry point via subcommands.

> The CLI application code is on the [`Cli_app`](https://github.com/caolei2-BGI/BRICK/tree/Cli_app) branch. Switch to that branch after cloning before installing or packaging.

## Command overview

| Subcommand | Description |
|--------|------|
| `qr` | Query entity relationships in the knowledge graph |
| `rk` | Rank or perform enrichment analysis on query results |
| `ga` | Read, write, convert, and merge graphs |
| `pl` | Static and interactive visualization |
| `emb` | Build an integrated UMAP graph from h5ad and CSV files |

```bash
brick --help              # List all subcommands
brick qr --help           # Show help for the qr subcommand
```

## Requirements

- Python >= 3.8 (pip installation)
- Network access (commands such as `qr` require access to the BRICK KG API)
- Optional: [Graphviz](https://graphviz.org/download/) (required to export SVG/PDF with `brick ga`)

---

## Option 1: Install with pip (recommended)

Suitable for users with an existing Python environment. After installation, the `brick` command is available globally.

### 1. Clone the repository and switch branches

The CLI code is on the **`Cli_app`** branch ([view branch](https://github.com/caolei2-BGI/BRICK/tree/Cli_app)):

```bash
# Option A: clone the specified branch directly (recommended)
git clone -b Cli_app https://github.com/caolei2-BGI/BRICK.git
cd BRICK

# Option B: clone first and then switch branches
git clone https://github.com/caolei2-BGI/BRICK.git
cd BRICK
git checkout Cli_app
```

### 2. Create a virtual environment (recommended)

```bash
# Linux / macOS
python3 -m venv .venv
source .venv/bin/activate

# Windows (PowerShell)
python -m venv .venv
.\.venv\Scripts\Activate.ps1
```

### 3. Install the package

**Editable install** (no reinstall required after code changes):

```bash
pip install --upgrade pip
pip install -e .
```

**Standard install**:

```bash
pip install .
```

### 4. Verify the installation

```bash
brick --help
brick qr --help
```

If the virtual environment's `Scripts` directory (Windows) or `bin` directory (Linux/macOS) has not been added to PATH, you can also run:

```bash
python cli.py --help
```

### 5. Optional isolated installation with pipx

To isolate the CLI from other Python packages, use [pipx](https://pipx.pypa.io/):

```bash
pipx install .
# Or from the cloned directory
pipx install -e /path/to/BRICK
```

### Usage examples

```bash
# Query entity relationships
brick qr "TP53" -o query.csv

# Ranking analysis
brick rk -i query.csv -o rank.csv

# Enrichment analysis
brick rk -i query.csv --enrich -o enrich.csv

# Read a graph file
brick ga read -i graph.graphml

# Interactive network visualization
brick pl interact -i graph.graphml -o network.html

# Build an integrated UMAP graph
brick emb -a data.h5ad -i query.csv -o output.h5ad
```

---

## Option 2: Build a binary (no Python environment required)

Suitable for distribution to users **without Python installed**. The build produces a standalone executable, still accessed through the `brick` command.

> **Notes**
> - The binary is usually large (typically 500 MB to 1 GB or more), because it includes scientific computing dependencies such as scanpy, scipy, and matplotlib
> - A binary built on Windows runs only on Windows; Linux and macOS require separate builds on their respective platforms
> - Exporting SVG/PDF still requires Graphviz on the user's machine

### 1. Clone the repository and switch branches

The CLI code is on the **`Cli_app`** branch ([view branch](https://github.com/caolei2-BGI/BRICK/tree/Cli_app)):

```bash
# Option A: clone the specified branch directly (recommended)
git clone -b Cli_app https://github.com/caolei2-BGI/BRICK.git
cd BRICK

# Option B: clone first and then switch branches
git clone https://github.com/caolei2-BGI/BRICK.git
cd BRICK
git checkout Cli_app
```

### 2. One-click build on Windows

Run the following in PowerShell:

```powershell
.\scripts\build_binary.ps1
```

The script automatically creates a `.venv-build` virtual environment, installs dependencies, and invokes PyInstaller.

### 3. Manual build (Windows / Linux / macOS)

```bash
# Create a dedicated build environment
python -m venv .venv-build

# Activate the environment
# Linux / macOS:
source .venv-build/bin/activate
# Windows PowerShell:
.\.venv-build\Scripts\Activate.ps1

# Install dependencies and PyInstaller
pip install --upgrade pip
pip install -e .
pip install pyinstaller

# Build the binary
pyinstaller brick.spec --noconfirm --clean
```

### 4. Build outputs

The default is **onedir** mode (faster startup and easier debugging):

```
dist/
└── brick/
    ├── brick          # Linux / macOS executable
    ├── brick.exe      # Windows executable
    └── _internal/     # Runtime dependencies (must be distributed with the executable)
```

Package the entire `dist/brick/` directory as a zip file. **Do not distribute only the exe file.**

### 5. Test the binary

```bash
# Linux / macOS
./dist/brick/brick --help
./dist/brick/brick qr --help

# Windows
.\dist\brick\brick.exe --help
.\dist\brick\brick.exe qr --help
```

### 6. Switch to single-file mode (optional)

To produce a single exe file, edit `brick.spec` and change `onedir` to `onefile` using the comments at the bottom of the file, then rerun `pyinstaller brick.spec --noconfirm --clean`.

Single-file mode starts more slowly because temporary files are unpacked at runtime, but it is easier to distribute.

---

## Dependencies

Core Python dependencies installed automatically by `pip install`:

- pandas, requests, scipy, networkx
- scanpy, matplotlib, pyvis, pydot

---

## FAQ

### The `brick` command is not found

Activate the virtual environment or reinstall in editable mode:

```bash
pip install -e .
```

### `brick ga` fails to export SVG/PDF

Install Graphviz and ensure `dot` is on the system PATH:

- Windows: download from [graphviz.org](https://graphviz.org/download/)
- Linux: `sudo apt install graphviz` or `sudo yum install graphviz`
- macOS: `brew install graphviz`

### `ModuleNotFoundError` after packaging

Add the missing modules to the `hiddenimports` list in `brick.spec` and rebuild.

### NumPy / pandas version conflicts during pip installation

Use a fresh virtual environment instead of mixing with an existing Anaconda environment:

```bash
python -m venv .venv
source .venv/bin/activate   # On Windows, run .venv\Scripts\Activate.ps1
pip install -e .
```

---

## License

See the LICENSE file in the repository, if present.
