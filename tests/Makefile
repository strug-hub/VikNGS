# Convenience wrapper for the VikNGS test suite.
# Run from the tests/ directory or: make -C tests <target>.

REPO := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))..)

.PHONY: test e2e stats unit build regen help

help:
	@echo "Targets:"
	@echo "  make build  — build bin/vikNGS, bin/sim-cli, bin/test-unit"
	@echo "  make test   — run all test layers (e2e + stats + unit)"
	@echo "  make e2e    — regression tests against golden files"
	@echo "  make stats  — simulation-driven statistical tests"
	@echo "  make unit   — Catch2 unit tests"
	@echo "  make regen  — regenerate e2e goldens (CONFIRM=1 required)"

build:
	$(MAKE) -C $(REPO)/bin
	$(MAKE) -C $(REPO)/bin sim-cli
	$(MAKE) -C $(REPO)/tests/unit

test: e2e stats unit

e2e:
	cd $(REPO) && python3 -m pytest tests/e2e -v

stats:
	cd $(REPO) && python3 -m pytest tests/stats -v

unit:
	$(MAKE) -C $(REPO)/tests/unit
	$(REPO)/bin/test-unit

regen:
ifeq ($(CONFIRM),1)
	cd $(REPO) && python3 -m tests.e2e.regen
else
	@echo "Refusing to overwrite goldens without CONFIRM=1."
	@echo "This will rewrite tests/e2e/golden/* from the current CLI output."
	@echo "Use:  make regen CONFIRM=1"
	@exit 1
endif
