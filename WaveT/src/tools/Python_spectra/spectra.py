#!/usr/bin/env python3
from config import read_input,validate
from dispatcher import CALCULATIONS

import sys
import time
from contextlib import redirect_stdout


def main():

    cfg = read_input(sys.argv[1])
    validate(cfg)
    print(f"Calculation : {cfg.calculation}")
    print(f"Medium      : {cfg.medium}")

    CALCULATIONS[cfg.calculation](cfg)


if __name__ == "__main__":

    start = time.time()

    with open(sys.argv[2], "w") as log:
        with redirect_stdout(log):
            main()
            print(f"Total time: {time.time()-start:.2f} s")
