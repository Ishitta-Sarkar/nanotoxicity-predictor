import sys
from pathlib import Path

SRC = Path(__file__).resolve().parent / "src"
if SRC.is_dir():
    sys.path.insert(0, str(SRC))

from nanotoxicity.analysis import main

if __name__ == "__main__":
    main()
