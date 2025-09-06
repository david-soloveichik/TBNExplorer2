#!/bin/bash

# Installation script for TBNExplorer2 shell completions

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COMPLETIONS_DIR="$SCRIPT_DIR/completions"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "TBNExplorer2 Shell Completions Installer"
echo "========================================="
echo

# Function to install zsh completions
install_zsh() {
    echo "Installing Zsh completions..."
    
    # Find the best location for zsh completions
    local install_dir=""
    local user_install=false
    
    # Check if we're in a virtual environment
    if [[ -n "$VIRTUAL_ENV" ]]; then
        install_dir="$VIRTUAL_ENV/share/zsh/site-functions"
        echo -e "${YELLOW}Detected virtual environment. Installing to: $install_dir${NC}"
    # Check for user-specific directory
    elif [[ -d "$HOME/.local/share/zsh/site-functions" ]]; then
        install_dir="$HOME/.local/share/zsh/site-functions"
        user_install=true
        echo -e "${YELLOW}Installing to user directory: $install_dir${NC}"
    elif [[ -d "$HOME/.zsh/completions" ]]; then
        install_dir="$HOME/.zsh/completions"
        user_install=true
        echo -e "${YELLOW}Installing to user directory: $install_dir${NC}"
    else
        # Create user directory if it doesn't exist
        install_dir="$HOME/.local/share/zsh/site-functions"
        user_install=true
        echo -e "${YELLOW}Creating and installing to: $install_dir${NC}"
    fi
    
    # Create directory if it doesn't exist
    mkdir -p "$install_dir"
    
    # Copy completion files
    cp "$COMPLETIONS_DIR/zsh/_tbnexplorer2" "$install_dir/"
    cp "$COMPLETIONS_DIR/zsh/_tbnexplorer2-filter" "$install_dir/"
    cp "$COMPLETIONS_DIR/zsh/_tbnexplorer2-ibot" "$install_dir/"
    
    echo -e "${GREEN}✓ Zsh completions installed successfully!${NC}"
    
    # Check if user needs to update their configuration
    if $user_install; then
        local zshrc="$HOME/.zshrc"
        local fpath_line="fpath=($install_dir \$fpath)"
        
        # Check if the directory is already in fpath in .zshrc
        if [[ -f "$zshrc" ]]; then
            if ! grep -q "$install_dir" "$zshrc"; then
                echo
                echo -e "${YELLOW}ACTION REQUIRED:${NC}"
                echo "To enable completions, add the following line to your ~/.zshrc"
                echo "(before the line that calls 'compinit'):"
                echo
                echo -e "${GREEN}    $fpath_line${NC}"
                echo
                echo "If you don't have compinit in your .zshrc, also add:"
                echo -e "${GREEN}    autoload -Uz compinit && compinit${NC}"
                echo
            else
                echo
                echo -e "${GREEN}✓ Directory already in fpath in ~/.zshrc${NC}"
            fi
        else
            echo
            echo -e "${YELLOW}ACTION REQUIRED:${NC}"
            echo "Create a ~/.zshrc file with these lines:"
            echo
            echo -e "${GREEN}    $fpath_line"
            echo "    autoload -Uz compinit && compinit${NC}"
            echo
        fi
    fi
    
    echo "To use the completions after updating .zshrc, run:"
    echo -e "${YELLOW}source ~/.zshrc${NC}"
    echo "or restart your terminal."
}

# Function to install for current shell session only
install_current_session() {
    echo "Installing completions for current session only..."
    
    if [[ -n "$ZSH_VERSION" ]]; then
        # Add to fpath for current session
        fpath=("$COMPLETIONS_DIR/zsh" $fpath)
        
        # Reload completions
        autoload -Uz compinit && compinit
        
        echo -e "${GREEN}✓ Completions loaded for current Zsh session!${NC}"
        echo "Note: This is temporary and will not persist after closing the terminal."
    else
        echo -e "${RED}Current shell is not Zsh. Please run this in a Zsh shell.${NC}"
        exit 1
    fi
}

# Parse command line arguments
case "${1:-}" in
    --current-session)
        install_current_session
        ;;
    --help|-h)
        echo "Usage: $0 [OPTIONS]"
        echo
        echo "Options:"
        echo "  --current-session    Install for current shell session only"
        echo "  --help, -h          Show this help message"
        echo
        echo "Without options, installs permanently to appropriate directory."
        ;;
    *)
        install_zsh
        ;;
esac