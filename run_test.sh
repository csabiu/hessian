#!/bin/bash
# Test script for Hessian density analysis program

echo "========================================="
echo "Hessian Density Analysis Test Suite"
echo "========================================="
echo ""

# Color codes for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Test parameters
TEST_FILE="test_input.dat"
BIN_SIZE="1.0"
CALC_DENSITY=".true."

# Check if test input file exists
if [ ! -f "$TEST_FILE" ]; then
    echo -e "${RED}ERROR: Test input file $TEST_FILE not found${NC}"
    exit 1
fi

echo -e "${YELLOW}Test Configuration:${NC}"
echo "  Input file: $TEST_FILE"
echo "  Bin size: $BIN_SIZE"
echo "  Calculate density: $CALC_DENSITY"
echo ""

# Count particles in input file
PARTICLE_COUNT=$(wc -l < "$TEST_FILE")
echo "  Particles in test file: $PARTICLE_COUNT"
echo ""

# Check if executable exists
if [ ! -f "./hessian" ]; then
    echo -e "${RED}ERROR: Hessian executable not found. Building...${NC}"
    ./make.sh
    if [ $? -ne 0 ]; then
        echo -e "${RED}ERROR: Build failed${NC}"
        exit 1
    fi
fi

echo -e "${YELLOW}Running Hessian analysis...${NC}"
echo ""

# Run the program
./hessian "$TEST_FILE" "$BIN_SIZE" "$CALC_DENSITY"

# Check if program completed successfully
if [ $? -eq 0 ]; then
    echo ""
    echo -e "${GREEN}=========================================${NC}"
    echo -e "${GREEN}Program completed successfully!${NC}"
    echo -e "${GREEN}=========================================${NC}"
    echo ""

    # Check for output files
    echo -e "${YELLOW}Checking output files:${NC}"

    OUTPUT_FILES=("density.dat" "smoothed_density.dat" "eigenvalue1.dat" "eigenvalue2.dat" "eigenvalue3.dat" "structure.dat")

    for file in "${OUTPUT_FILES[@]}"; do
        if [ -f "$file" ]; then
            SIZE=$(stat -f%z "$file" 2>/dev/null || stat -c%s "$file" 2>/dev/null)
            echo -e "  ${GREEN}✓${NC} $file (${SIZE} bytes)"
        else
            echo -e "  ${YELLOW}⊗${NC} $file (not created - may be expected)"
        fi
    done

    echo ""
    echo -e "${GREEN}Test PASSED${NC}"
    exit 0
else
    echo ""
    echo -e "${RED}=========================================${NC}"
    echo -e "${RED}Program failed with error code $?${NC}"
    echo -e "${RED}=========================================${NC}"
    echo ""
    echo -e "${RED}Test FAILED${NC}"
    exit 1
fi
