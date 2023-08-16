#!/usr/bin/env python3

import argparse
import pandas as pd
from pathlib import Path
import sqlite3

from supermodel.get_dbs import get_dbs, apply


def main(tsv_dir, sql_file):
    # Put all the tables in one sqlite database. First remove the previous db
    # if it exists
    Path(sql_file).unlink(missing_ok=True)
    con = sqlite3.connect(sql_file)
    cur = con.cursor()

    for name, data in get_dbs(tsv_dir).items():
        if name == "bigg_db_network_m":
            data = data.pipe(apply, "reactions", lambda x: ",".join(x))
        data.to_sql(name, con)

    con.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--tsv_dir", help="Path to all tsv files (default=data)", default="data"
    )
    parser.add_argument(
        "--sql_file",
        help="Path to output sqlite file (default=data/all.sqlite",
        default="data/all.sqlite",
    )
    args = parser.parse_args()

    main(**vars(args))
