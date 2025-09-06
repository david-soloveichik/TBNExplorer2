# Shell Completions for TBNExplorer2

This project includes zsh shell completions for all TBNExplorer2 CLI tools.

## Features

The completion scripts provide intelligent tab completion for:

- **Command options and flags**: All command-line options are available via tab completion
- **File paths**: Automatic completion for relevant file types:
  - `.tbn` files for input
  - `.tbnpolys` files for polymer specifications
  - `.tbnpolymat` files (automatically inferred)
  - `.txt` files for constraints
- **Value completions**: 
  - Concentration units (nM, pM, uM, mM, M)
  - Parameter formats for `--parametrized` option
- **Dynamic completions**:
  - Monomer names extracted from `.tbn` files (for `tbnexplorer2-filter`)

## Installation

### Quick Install (Recommended)

Run the installation script:

```bash
./install-completions.sh
```

This will:
1. Install completions to the appropriate directory
2. Update your `~/.zshrc` if needed
3. Provide instructions to activate completions

To activate immediately:
```bash
source ~/.zshrc
```

### Manual Installation

1. Copy the completion files to your zsh completions directory:
   ```bash
   cp completions/zsh/_tbnexplorer2* ~/.local/share/zsh/site-functions/
   ```

2. Add the directory to your `fpath` in `~/.zshrc`:
   ```bash
   fpath=(~/.local/share/zsh/site-functions $fpath)
   ```

3. Ensure `compinit` is called in your `~/.zshrc`:
   ```bash
   autoload -Uz compinit && compinit
   ```

4. Reload your shell configuration:
   ```bash
   source ~/.zshrc
   ```

### Installation with pip/setup.py

When installing the package with pip, completions are automatically installed to the system directory:

```bash
pip install .
```

## Usage Examples

### tbnexplorer2
```bash
# Tab completion for input files
tbnexplorer2 <TAB>           # Shows .tbn files

# Tab completion for options
tbnexplorer2 --<TAB>         # Shows all available options

# Tab completion for parametrized variables
tbnexplorer2 file.tbn --parametrized <TAB>   # Shows format hint
```

### tbnexplorer2-filter
```bash
# Tab completion for monomer names (extracted from .tbn file)
tbnexplorer2-filter example.tbn <TAB>    # Shows available monomer names

# Tab completion for constraints file
tbnexplorer2-filter example.tbn --constraints-file <TAB>   # Shows .txt files
```

### tbnexplorer2-ibot
```bash
# Tab completion for on-target polymers file
tbnexplorer2-ibot input.tbn <TAB>        # Shows .tbnpolys files

# Tab completion for units
tbnexplorer2-ibot input.tbn target.tbnpolys --generate-tbn 100 <TAB>  # Shows unit options
```

## Troubleshooting

### Completions not working

1. Verify completions are installed:
   ```bash
   ls ~/.local/share/zsh/site-functions/_tbnexplorer2*
   ```

2. Check that the directory is in your fpath:
   ```bash
   echo $fpath | grep -o '[^ ]*zsh/site-functions'
   ```

3. Rebuild completion cache:
   ```bash
   rm -f ~/.zcompdump && compinit
   ```

4. Ensure you're using zsh:
   ```bash
   echo $SHELL
   ```

### Updating completions

After updating the completion scripts, rebuild the completion cache:
```bash
rm -f ~/.zcompdump && compinit
```

## Development

The completion scripts are located in `completions/zsh/`:
- `_tbnexplorer2`: Main tool completions
- `_tbnexplorer2-filter`: Filter tool completions  
- `_tbnexplorer2-ibot`: IBOT algorithm tool completions

To test changes:
1. Edit the completion script
2. Copy to your completions directory
3. Reload completions: `unfunction _tbnexplorer2 2>/dev/null; autoload -Uz _tbnexplorer2`
4. Test with tab completion

## Contributing

When adding new command-line options to the tools, please update the corresponding completion script in `completions/zsh/`.