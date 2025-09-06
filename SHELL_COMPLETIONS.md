# Shell Completions for TBNExplorer2

Tab completion for all TBNExplorer2 commands, automatically generated from the command definitions.

## Quick Start

**Choose your installation method:**

- **Option A: Dynamic Loading** (Recommended for most users)  
  Add a few lines to your shell config file. 

- **Option B: Static Files** (For zsh users who prefer not editing config files)  
  Install completion files once. 

Both methods provide the same features and automatically reflect any changes to command options.

## Features

- Complete command options (`--output`, `--verbose`, etc.)
- Smart file suggestions (`.tbn`, `.tbnpolys`, `.txt` files)
- Concentration units (`nM`, `pM`, `uM`, `mM`, `M`)
- Dynamic monomer name extraction from `.tbn` files
- Parameter format hints for `--parametrized`

## Requirements

- Python 3.8+
- argcomplete package (installed automatically with TBNExplorer2)
- bash 4.2+ or zsh

**macOS Users:** The default bash is too old (3.2). Use zsh instead (default on macOS 10.15+) or install newer bash with `brew install bash`.

## Installation

### Option A: Dynamic Loading (Recommended)

Add these lines to your shell configuration file, then reload:

**For bash** (`~/.bashrc`):
```bash
# TBNExplorer2 completions
eval "$(register-python-argcomplete tbnexplorer2)"
eval "$(register-python-argcomplete tbnexplorer2-filter)"
eval "$(register-python-argcomplete tbnexplorer2-ibot)"
```

**For zsh** (`~/.zshrc`):
```zsh
# TBNExplorer2 completions
autoload -U bashcompinit && bashcompinit
eval "$(register-python-argcomplete tbnexplorer2)"
eval "$(register-python-argcomplete tbnexplorer2-filter)"
eval "$(register-python-argcomplete tbnexplorer2-ibot)"
```

Then reload your shell:
```bash
source ~/.bashrc  # or ~/.zshrc for zsh
```

### Option B: Static Files (zsh only)

**ONE-TIME INSTALLATION** - Run this once and forget about it:

```bash
# Install static completion files
./activate-completions.sh --install-zsh-static

# Reload zsh
exec zsh
```

This creates completion files in `~/.local/share/zsh/site-functions/`. These files are smart - they call the Python scripts to generate completions, so they automatically stay current with any code changes.

**Important:** You do NOT need to rerun this command when parameters change or when you update the code. The completions will automatically reflect all changes.

### Helper Script

The `activate-completions.sh` script can help with setup:

```bash
# Show instructions for dynamic loading
./activate-completions.sh --per-command

# Install static zsh files (one-time only)
./activate-completions.sh --install-zsh-static

# Test if completions are working
./activate-completions.sh --test
```

## Usage Examples

Once installed, use tab completion with any command:

```bash
# Complete command options
tbnexplorer2 --<TAB>                     # Shows all options

# Complete file types
tbnexplorer2 <TAB>                       # Shows .tbn files
tbnexplorer2 file.tbn --output <TAB>     # Shows .tbnpolys files

# Complete monomer names (extracted from the .tbn file)
tbnexplorer2-filter example.tbn <TAB>    # Shows monomer names

# Complete concentration units
tbnexplorer2-ibot input.tbn target.tbnpolys --generate-tbn 100 <TAB>
```

## Testing Your Setup

Quick test:
```bash
tbnexplorer2 --<TAB><TAB>
```

You should see a list of available options. If not, see Troubleshooting below.

## Troubleshooting

### Completions not working?

1. **Check argcomplete is installed:**
   ```bash
   python3 -c "import argcomplete" && echo "OK"
   ```
   If this fails: `pip install argcomplete`

2. **For dynamic loading:** Verify the eval commands are in your shell config and you've reloaded
   
3. **For static files:** Check that `~/.local/share/zsh/site-functions/` exists and contains `_tbnexplorer2*` files

4. **Reload your shell:**
   ```bash
   exec $SHELL  # Start fresh shell
   ```

### macOS Issues

If completions don't work on macOS:
```bash
# Check your bash version
echo $BASH_VERSION

# If it shows 3.2, switch to zsh:
chsh -s /bin/zsh
```

## For Developers

### Adding New Commands

Completions are automatically generated from argparse definitions. Just ensure:

1. Your script has the marker after the shebang:
   ```python
   #!/usr/bin/env python3
   # PYTHON_ARGCOMPLETE_OK
   ```

2. Enable argcomplete:
   ```python
   import argcomplete
   # ... create parser ...
   argcomplete.autocomplete(parser)
   ```

### Custom Completers

For specialized completions, add a completer to `tbnexplorer2/completers.py`:

```python
from .completers import MyCustomCompleter

argument = parser.add_argument("--my-option")
argument.completer = MyCustomCompleter
```

### Important Notes

- **All completion methods auto-update:** Whether using dynamic loading or static files, completions automatically reflect code changes
- **Static files are wrappers:** They don't contain the actual completions, just calls to generate them
- **No manual regeneration needed:** Changes to command options are picked up automatically

## Summary

1. Install completions once using either method
2. Completions automatically stay current with code changes
3. Use `<TAB>` to complete commands, files, and options
4. No maintenance required after initial setup