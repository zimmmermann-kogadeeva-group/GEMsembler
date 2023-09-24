#!/usr/bin/env python3

import argparse
import pandas as pd
from pathlib import Path
import sqlite3

from supermodel.dbs import (
    apply,
    seed_orig_m,
    seed_orig_r,
    seed_addit_m,
    seed_addit_r,
    kegg_bigg_m,
    kegg_bigg_r,
    bigg_all_m,
    bigg_all_r,
    bigg_db_network_m,
    bigg_db_network_r,
)


def main(sql_file):
    # Put all the tables in one sqlite database. First remove the previous db
    # if it exists
    Path(sql_file).unlink(missing_ok=True)
    with sqlite3.connect(sql_file) as con:
        seed_orig_m.to_sql("seed_orig_m", con)
        seed_orig_r.to_sql("seed_orig_r", con)
        seed_addit_m.to_sql("seed_addit_m", con)
        seed_addit_r.to_sql("seed_addit_r", con)
        kegg_bigg_m.to_sql("kegg_bigg_m", con)
        kegg_bigg_r.to_sql("kegg_bigg_r", con)
        bigg_all_m.to_sql("bigg_all_m", con)
        bigg_all_r.to_sql("bigg_all_r", con)
        bigg_db_network_m.pipe(apply, "reactions", lambda x: ",".join(x)).to_sql(
            "bigg_db_network_m", con
        )
        bigg_db_network_r.to_sql("bigg_db_network_r", con)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--sql_file",
        help="Path to output sqlite file (default=data/all.sqlite",
        default="data/all.sqlite",
    )
    args = parser.parse_args()

    main(**vars(args))
