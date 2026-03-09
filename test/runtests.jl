using TestItemRunner

# Standard package test entrypoint. Domain-specific thin wrappers live in this
# directory and forward to tests/<family>/runtests.jl.
@run_package_tests