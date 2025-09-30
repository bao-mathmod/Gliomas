#!/bin/bash

# Define the list of tools to check
tools=(
    curl
    column
    awk
    grep
    sed
    wc
    prefetch
    fasterq-dump
    pigz
    zcat
)

# Initialize a counter for missing tools
missing_tools_count=0

# Loop through each tool in the list
for tool in "${tools[@]}"; do
    # Use 'command -v' to check if the tool is in the system's PATH
    if command -v "$tool" &> /dev/null; then
        echo "✅ $tool is installed."
    else
        echo "❌ $tool is NOT installed."
        # Increment the counter for missing tools
        ((missing_tools_count++))
    fi
done

echo "" # Add a newline for better readability

# Provide a summary of the results
if [ "$missing_tools_count" -eq 0 ]; then
    echo "🎉 All specified tools are installed. Your environment is ready!"
else
    echo "⚠️  $missing_tools_count tool(s) are missing. Please install them to proceed."
fi