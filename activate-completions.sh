#!/bin/bash

# Activation script for TBNExplorer2 argcomplete-based shell completions

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo "TBNExplorer2 Shell Completions Activation"
echo "=========================================="
echo

# Check if argcomplete is installed
check_argcomplete() {
    if python3 -c "import argcomplete" 2>/dev/null; then
        return 0
    else
        return 1
    fi
}

# Function to get user's shell type (prefer login shell from $SHELL)
get_shell() {
    if [[ -n "$SHELL" ]]; then
        basename "$SHELL"
        return
    fi

    # Fallbacks if $SHELL is not set
    if command -v ps >/dev/null 2>&1; then
        ps -p $$ -o comm= 2>/dev/null | awk -F/ '{print $NF}'
        return
    fi

    if [[ -n "$ZSH_VERSION" ]]; then
        echo "zsh"
        return
    fi
    if [[ -n "$BASH_VERSION" ]]; then
        echo "bash"
        return
    fi

    echo "unknown"
}

# Main installation logic
main() {
    local shell_type=$(get_shell)
    local method="${1:-per-command}"
    
    echo -e "${BLUE}Detected shell: $shell_type${NC}"
    echo
    
    # Check argcomplete installation
    if ! check_argcomplete; then
        echo -e "${RED}Error: argcomplete is not installed${NC}"
        echo
        echo "Please install argcomplete first:"
        echo -e "${GREEN}pip install argcomplete${NC}"
        echo
        echo "Or if you've installed tbnexplorer2:"
        echo -e "${GREEN}pip install -e .${NC}"
        exit 1
    fi
    
    echo -e "${GREEN}✓ argcomplete is installed${NC}"
    echo
    
    case "$method" in
        per-command|--per-command)
            echo "Setting up per-command argcomplete registration..."
            echo
            
            local rc_file=""
            case "$shell_type" in
                bash)
                    rc_file="$HOME/.bashrc"
                    ;;
                zsh)
                    rc_file="$HOME/.zshrc"
                    ;;
                *)
                    echo -e "${RED}Unsupported shell type: $shell_type${NC}"
                    exit 1
                    ;;
            esac
            
            echo "Add the following lines to your $rc_file:"
            echo
            if [[ "$shell_type" == "zsh" ]]; then
                cat <<'EOF'
# TBNExplorer2 completions (zsh)
autoload -U bashcompinit
bashcompinit
eval "$(register-python-argcomplete tbnexplorer2)"
eval "$(register-python-argcomplete tbnexplorer2-filter)"
eval "$(register-python-argcomplete tbnexplorer2-ibot)"
EOF
            else
                cat <<'EOF'
# TBNExplorer2 completions (bash)
eval "$(register-python-argcomplete tbnexplorer2)"
eval "$(register-python-argcomplete tbnexplorer2-filter)"
eval "$(register-python-argcomplete tbnexplorer2-ibot)"
EOF
            fi
            echo
            echo "Then reload your shell configuration:"
            echo -e "${YELLOW}source $rc_file${NC}"
            ;;
            
        install-zsh-static|--install-zsh-static)
            echo "Installing zsh completion wrapper files..."
            echo
            
            # Ensure register tool is available
            if ! command -v register-python-argcomplete >/dev/null 2>&1; then
                echo -e "${RED}Error: register-python-argcomplete not found${NC}"
                echo "Install argcomplete (pip install argcomplete) and ensure your Python environment is active."
                exit 1
            fi

            # Target directory for zsh completion functions
            local target_dir="$HOME/.local/share/zsh/site-functions"
            mkdir -p "$target_dir"

            local scripts=("tbnexplorer2" "tbnexplorer2-filter" "tbnexplorer2-ibot")
            local any_installed=false

            for script in "${scripts[@]}"; do
                if command -v "$script" >/dev/null 2>&1; then
                    local outfile="$target_dir/_$script"
                    if register-python-argcomplete --shell zsh "$script" > "$outfile" 2>/dev/null; then
                        echo -e "${GREEN}✓ Installed: $outfile${NC}"
                        any_installed=true
                    else
                        echo -e "${YELLOW}✗ Failed to generate completion for $script${NC}"
                    fi
                else
                    echo -e "${YELLOW}? $script: not found in PATH (skipping)${NC}"
                fi
            done

            echo
            if $any_installed; then
                echo "If this directory is in your fpath before compinit, completions will be active:"
                echo "  $target_dir"
                echo
                echo "Reload zsh completions:"
                echo "  exec zsh"
                echo "  # or: rm -f ~/.zcompdump-$HOST && compinit"
            else
                echo -e "${YELLOW}No completion files were installed. Ensure commands are installed and in PATH.${NC}"
                exit 1
            fi
            ;;
            
        test|--test)
            echo "Testing argcomplete setup (per-command)..."
            echo
            
            local scripts=("tbnexplorer2" "tbnexplorer2-filter" "tbnexplorer2-ibot")
            local all_ok=true

            # Check register tool
            if command -v register-python-argcomplete >/dev/null 2>&1; then
                echo -e "${GREEN}✓ register-python-argcomplete is available${NC}"
            else
                echo -e "${RED}✗ register-python-argcomplete not found${NC}"
                all_ok=false
            fi

            # Probe per-command registration generation (dry-run to stdout)
            for script in "${scripts[@]}"; do
                if command -v "$script" >/dev/null 2>&1; then
                    if register-python-argcomplete "$script" >/dev/null 2>&1; then
                        echo -e "${GREEN}✓ $script: registration probe OK${NC}"
                    else
                        echo -e "${YELLOW}✗ $script: registration probe failed${NC}"
                        all_ok=false
                    fi
                else
                    echo -e "${YELLOW}? $script: not found in PATH${NC}"
                fi
            done

            echo
            if $all_ok; then
                echo -e "${GREEN}Per-command setup looks good.${NC}"
                echo "Add eval lines to your shell RC (see --per-command), then test:"
                echo "  tbnexplorer2 --<TAB><TAB>"
            else
                echo -e "${YELLOW}Some checks failed.${NC}"
                echo "Ensure argcomplete is installed and try per-command setup (see --per-command)."
            fi
            ;;
            
        help|--help|-h)
            echo "Usage: $0 [METHOD]"
            echo
            echo "Methods:"
            echo "  per-command, --per-command   Set up per-command registration (recommended)"
            echo "  install-zsh-static, --install-zsh-static  Install zsh completion files under ~/.local/share/zsh/site-functions"
            echo "  test, --test            Test if argcomplete is properly configured"
            echo "  help, --help, -h        Show this help message"
            echo
            echo "Examples:"
            echo "  $0 --per-command       # Per-command registration instructions (recommended)"
            echo "  $0 --install-zsh-static  # Install zsh completion wrapper files"
            echo "  $0 --test              # Test the setup"
            echo
            echo "Requirements:"
            echo "  - Python 3.8+"
            echo "  - argcomplete package (pip install argcomplete)"
            echo "  - bash 4.2+ or zsh"
            ;;
            
        *)
            echo -e "${RED}Unknown method: $method${NC}"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
}

# Parse arguments and run
main "$@"
