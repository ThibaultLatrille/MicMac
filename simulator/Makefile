
# ==============================================================================================================
#  COMPILATION
# ==============================================================================================================
.PHONY: all # Requires: cmake 3.1.0 or better
all: build
	@cd build ; make --no-print-directory -j8

build: CMakeLists.txt
	@rm -rf build
	@mkdir build
	@cd build ; cmake ..

.PHONY: coverage
coverage:
	@rm -rf build
	@mkdir build
	@cd build ; cmake -DCOVERAGE_MODE=ON ..
	@make --no-print-directory test

.PHONY: debug
debug:
	@rm -rf build
	@mkdir build
	@cd build ; cmake -DDEBUG_MODE=ON ..
	@make --no-print-directory

.PHONY: release
release:
	@rm -rf build
	@mkdir build
	@cd build ; cmake ..
	@make --no-print-directory

.PHONY: clean
clean:
	@rm -rf build
	@rm -rf test

# ==============================================================================================================
#  CODE QUALITY
# ==============================================================================================================
.PHONY: format # Requires: clang-format
format:
	@clang-format -i `find src/ -name *.*pp`

# ==============================================================================================================
#  TESTING
# ==============================================================================================================

.PHONY: test
test: debug
	@rm -rf test
	@mkdir test
	@echo "\n\e[35m\e[1m== Neutral run =========================================================\e[0m"
	build/neutral --output test/neutral_run
	@echo "\n\e[35m\e[1m== Stabilizing run =====================================================\e[0m"
	build/stabilizing --output test/stabilizing_run
	@echo "\n\e[35m\e[1m== Directional run =====================================================\e[0m"
	build/directional --output test/directional_run
	@echo "\n\e[35m\e[1m== Moving optimum run =====================================================\e[0m"
	build/moving_optimum --output test/moving_optimum_run