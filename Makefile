#==============================================================================
# IOSEE - Infrared Ocean Surface Effective Emissivity
# Top-Level Makefile with Parallel Build Support
#==============================================================================

# Auto-detect CPU cores for parallel builds
NPROC := $(shell nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

# Auto-detect version from VERSION.md (table row: "| Version | **X.Y.Z** |")
VERSION := $(shell grep -E '^\| Version' VERSION.md 2>/dev/null | sed -E 's/.*\*\*([0-9]+\.[0-9]+\.[0-9]+)\*\*.*/\1/' || echo "unknown")

# Build directory
BUILD_DIR := build

# CMake generator (use Ninja if available for faster builds)
CMAKE_GENERATOR := $(shell command -v ninja >/dev/null 2>&1 && echo "Ninja" || echo "Unix Makefiles")

.PHONY: release debug clean reproducible install help openmp-info

# Default target
all: release

# Release build with parallel compilation
release:
	@echo "Building Ocean Emissivity (Release mode) with $(NPROC) cores..."
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && cmake -G "$(CMAKE_GENERATOR)" -DCMAKE_BUILD_TYPE=Release ..
	@cmake --build $(BUILD_DIR) -j$(NPROC)
	@echo "                   "
	@echo "The Infrared Ocean Surface Effective Emissivity (IOSEE) package version $(VERSION) build successfully!"
	@echo "                   "
	@echo "========================================================="
	@echo "Copyright (C) 2025 - :"
	@echo "     Texas A&M University"
	@echo "      "
	@echo "Maintainer (contact):  "
	@echo "     Atmospheric & Oceanic Optics Group, Department of Atmospheric Sciences, Texas A&M University"
	@echo "      "
	@echo "Citation:          "
	@echo "     For any application using the IOSEE package, please cite Wei et al. (2026)."
	@echo "     Wei J., Yang, P. (2026). IOSEE: An Efficient and Flexible Infrared Ocean"
	@echo "                              Surface Effective Emissivity Package for Remote Sensing,"
	@echo "                              Numerical Weather Prediction, and Climate Modeling Applications."
	@echo "                              Journal of Quantitative Spectroscopy and Radiative Transfer. ***, ****. "	
	@echo "========================================================="
	@echo "                   "
	@echo "========================================================="
	@echo "How to run"
	@echo "                   "
	@echo "The folder 'data' in the ./examples/ is a symlink that points to ./data"
	@echo "The command 'iosee' in the ./examples/ is a symlink that points to ./build/iosee"
	@echo "   "
	@echo "Executable command (in the ./examples):"
	@echo "    cd examples"
	@echo "    ./iosee user_defined.config output=user_defined.nc"	
	@echo "                         "
	@echo "========================================================="

# Debug build with parallel compilation
debug:
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && cmake -G "$(CMAKE_GENERATOR)" -DCMAKE_BUILD_TYPE=Debug ..
	@cmake --build $(BUILD_DIR) -j$(NPROC)
	@echo ""
	@echo "IOSEE version $(VERSION) debug build complete. Executable: $(BUILD_DIR)/iosee"

# Reproducible build with strict IEEE 754 compliance
reproducible:
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && cmake -G "$(CMAKE_GENERATOR)" -DCMAKE_BUILD_TYPE=Release -DENABLE_REPRODUCIBLE=ON ..
	@cmake --build $(BUILD_DIR) -j$(NPROC)
	@echo ""
	@echo "IOSEE version $(VERSION) reproducible build complete (strict IEEE 754). Executable: $(BUILD_DIR)/iosee"

# Build without OpenMP (serial mode)
serial:
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && cmake -G "$(CMAKE_GENERATOR)" -DCMAKE_BUILD_TYPE=Release -DENABLE_OPENMP=OFF ..
	@cmake --build $(BUILD_DIR) -j$(NPROC)
	@echo ""
	@echo "IOSEE version $(VERSION) serial build complete. Executable: $(BUILD_DIR)/iosee"

# Clean build artifacts
clean:
	@echo "Cleaning build files..."
	# Remove the build directories
	@rm -rf $(BUILD_DIR)
	# Remove the symbolic links to the executable and data
	@rm -f examples/iosee examples/data
	# Remove any stray object files
	@find . -name "*.o" -type f -delete
	# Remove any module files
	@find . -name "*.mod" -type f -delete
	# Remove any output files
	@find . -name "*.out" -type f -delete
	# Remove any log files
	@find . -name "*.log" -type f -delete
	@echo "Clean complete"

# Full clean including CMake cache
distclean: clean
	@rm -f CMakeCache.txt cmake_install.cmake
	@rm -rf CMakeFiles
	@echo "Distribution clean complete"

# Install to system (requires sudo on most systems)
install: release
	@cmake --install $(BUILD_DIR)

# Show build information
info:
	@echo "IOSEE Build System"
	@echo "==================="
	@echo "Version:            $(VERSION)"
	@echo "Detected CPU cores: $(NPROC)"
	@echo "CMake generator:    $(CMAKE_GENERATOR)"
	@echo "Build directory:    $(BUILD_DIR)"
	@echo ""
	@echo "Targets:"
	@echo "  make release      - Optimized release build (default)"
	@echo "  make debug        - Debug build with symbols"
	@echo "  make reproducible - Bit-for-bit reproducible build"
	@echo "  make serial       - Build without OpenMP"
	@echo "  make clean        - Remove build artifacts"
	@echo "  make install      - Install to system"

# Show OpenMP configuration hints
openmp-info:
	@echo "OpenMP Runtime Configuration"
	@echo "============================"
	@echo "Recommended settings for optimal performance:"
	@echo ""
	@echo "  export OMP_NUM_THREADS=$(NPROC)"
	@echo "  export OMP_PROC_BIND=spread"
	@echo "  export OMP_PLACES=cores"
	@echo ""
	@echo "Or use the provided runtime wrapper:"
	@echo "  ./run_iosee.sh examples/spectral_option1.config"

help: info
