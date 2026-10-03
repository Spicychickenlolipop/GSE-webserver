"""
Quick standalone check: is R + rpy2 + Bioconductor limma set up correctly?

Usage:
    python check_r_setup.py

Doesn't touch the job queue or GEO — just verifies the R bridge works, so you
can debug R/rpy2 installation issues independently of the rest of the app.
"""

import sys

import r_bridge

available, message = r_bridge.check_available()

if available:
    print(f"✓ {message}")
    sys.exit(0)
else:
    print(f"✗ R/limma is NOT available: {message}", file=sys.stderr)
    print(
        "\nThe pipeline will still run using the Python fallback (quantile "
        "normalization + Welch's t-test) — this just means microarray "
        "results won't use the field-standard limma statistics.\n"
        "\nTo enable limma:\n"
        "  1. Install R: https://cran.r-project.org/\n"
        "  2. In R: install.packages('BiocManager'); BiocManager::install('limma')\n"
        "  3. pip install rpy2\n",
        file=sys.stderr,
    )
    sys.exit(1)
