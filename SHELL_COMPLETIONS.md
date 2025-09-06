# Shell Completions for TBNExplorer2

This project uses [argcomplete](https://github.com/kislyuk/argcomplete) to provide intelligent tab completion for all TBNExplorer2 CLI tools.

## Features

The argcomplete-based completion system provides:

- Automatic completion generation from argparse definitions
- Multi-shell support: bash (4.2+) and zsh
- Intelligent completions:
  - Command options and flags with descriptions
  - File path completion for specific file types (`.tbn`, `.tbnpolys`, `.txt`)
  - Concentration units (nM, pM, uM, mM, M)
  - Dynamic monomer name extraction from `.tbn` files
  - Parameter format hints for `--parametrized` option
- Always up-to-date: completions reflect any changes to command-line arguments

## Requirements

- Python 3.8+
- argcomplete package (automatically installed with tbnexplorer2)
- bash 4.2+ or zsh

### macOS Users
macOS ships with an outdated version of bash (3.2). You can either:
1. Use zsh (recommended, default on macOS 10.15+)
2. Install a newer bash version: `brew install bash`

## Installation

### Method 1: Per-Command Registration (Recommended)

Register each command in your shell configuration:

For bash, add to `~/.bashrc`:

```bash
eval "$(register-python-argcomplete tbnexplorer2)"
eval "$(register-python-argcomplete tbnexplorer2-filter)"
eval "$(register-python-argcomplete tbnexplorer2-ibot)"

# Reload your shell configuration
source ~/.bashrc
```

For zsh, add to `~/.zshrc` (bash completion emulation is required):

```zsh
autoload -U bashcompinit
bashcompinit
eval "$(register-python-argcomplete tbnexplorer2)"
eval "$(register-python-argcomplete tbnexplorer2-filter)"
eval "$(register-python-argcomplete tbnexplorer2-ibot)"

# Reload your shell configuration
source ~/.zshrc
```

### Method 2: Static zsh Completion Files (No RC edits)

If your `~/.zshrc` already adds a directory to `fpath` and runs `compinit`, you can drop static completion files there. A common user path is `~/.local/share/zsh/site-functions`.

```bash
mkdir -p ~/.local/share/zsh/site-functions
register-python-argcomplete --shell zsh tbnexplorer2 > ~/.local/share/zsh/site-functions/_tbnexplorer2
register-python-argcomplete --shell zsh tbnexplorer2-filter > ~/.local/share/zsh/site-functions/_tbnexplorer2-filter
register-python-argcomplete --shell zsh tbnexplorer2-ibot > ~/.local/share/zsh/site-functions/_tbnexplorer2-ibot

# Reload zsh completions
exec zsh   # or: rm -f ~/.zcompdump-$(hostname) && compinit
```

Ensure the target directory is in your `fpath` before `compinit` runs.

Alternatively, use the helper script to install these files automatically:

```bash
./activate-completions.sh --install-zsh-static
```

### Method 3: Using pip Installation

When installing with pip, argcomplete is automatically installed as a dependency:

```bash
pip install .
# or for development
pip install -e .
```

Then activate completions using one of the methods above.

## Usage

Once activated, use tab completion with any TBNExplorer2 command:

### tbnexplorer2
```bash
# Complete input files
tbnexplorer2 <TAB>                     # Shows .tbn files

# Complete options
tbnexplorer2 --<TAB>                   # Shows all available options

# Complete file outputs
tbnexplorer2 file.tbn --output <TAB>   # Shows .tbnpolys files

# Complete parametrized variables
tbnexplorer2 file.tbn --parametrized <TAB>  # Shows format hints
```

### tbnexplorer2-filter
```bash
# Complete monomer names (dynamically extracted from .tbn file)
tbnexplorer2-filter example.tbn <TAB>  # Shows available monomer names

# Complete constraints file
tbnexplorer2-filter example.tbn --constraints-file <TAB>  # Shows .txt files
```

### tbnexplorer2-ibot
```bash
# Complete on-target polymers file
tbnexplorer2-ibot input.tbn <TAB>      # Shows .tbnpolys files

# Complete concentration units
tbnexplorer2-ibot input.tbn target.tbnpolys --generate-tbn 100 <TAB>  # Shows units
```

## Testing Completions

To verify completions are set up:

```bash
# Run the test command (probes per-command registration)
./activate-completions.sh --test

# Manual interactive test
tbnexplorer2 --<TAB><TAB>  # Should show available options
```

## Troubleshooting

### Completions Not Working

1. **Check argcomplete is installed**:
   ```bash
   python3 -c "import argcomplete"
   ```

2. **Verify activation**:
   - For per-command: Ensure `register-python-argcomplete` eval commands are in your shell RC file (and `bashcompinit` for zsh)
   - For static zsh files: Ensure the directory is in `fpath` before `compinit`

3. **Check shell compatibility**:
   ```bash
   echo $BASH_VERSION  # Should be 4.2+
   echo $ZSH_VERSION   # Any version
   ```

4. **Reload shell configuration**:
   ```bash
   source ~/.bashrc  # or ~/.zshrc
   # Or start a new terminal
   ```

### macOS-Specific Issues

If using the system bash (3.2):
```bash
# Switch to zsh
chsh -s /bin/zsh

# Or install newer bash
brew install bash
echo '/opt/homebrew/bin/bash' | sudo tee -a /etc/shells
chsh -s /opt/homebrew/bin/bash
```

## How It Works

Argcomplete works by:
1. Intercepting tab completion requests from the shell
2. Running the Python script in completion mode to generate suggestions
3. Returning appropriate completions based on context

Per-command registration and static zsh files do not require wrapper markers.

The custom completers in `tbnexplorer2/completers.py` provide specialized completion logic for:
- File types with specific extensions
- Dynamic extraction of monomer names from .tbn files
- Concentration units and parameter formats

## Development

When adding new command-line options:

1. **No manual updates needed**: Completions are automatically generated from argparse
2. **For custom completions**: Add a completer function to `tbnexplorer2/completers.py`
3. **Assign the completer**: Set `argument.completer = your_completer_function`

Example:
```python
# In your CLI file
from .completers import MyCustomCompleter

arg = parser.add_argument("--my-option")
arg.completer = MyCustomCompleter
```

## Contributing

When adding new CLI tools or modifying arguments:
- Ensure the `# PYTHON_ARGCOMPLETE_OK` marker is present after the shebang
- Import and enable argcomplete: `argcomplete.autocomplete(parser)`
- Add custom completers as needed in `tbnexplorer2/completers.py`
- Test completions work with the new arguments
