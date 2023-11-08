import json
from functools import wraps
import pandas as pd
from pathlib import Path
import re


# helper functions for pandas dataframes
def separate(data, col, into=None, sep=" ", **kwargs):
    into = into or col
    data[into] = data[col].str.split(sep, **kwargs)
    return data


def replace(data, col, pat, repl, regex=True, **kwargs):
    data[col] = data[col].str.replace(pat, repl, regex=regex)
    return data


def set_value(data, col, value):
    data[col] = value
    return data


def colapply(data, col, func):
    data[col] = func(data[col])
    return data


def apply(data, col, func, new_col=None, **kwargs):
    new_col = new_col or col
    data[new_col] = data[col].apply(func, **kwargs)
    return data


def cache_file(func):
    """
    Decorator for caching mapping dictionaries in ~/.gemsembler
    """

    @wraps(func)
    def wrapper_decorator(*args, **kwargs):
        # Create the directory holding
        cache_path = Path("~").expanduser() / ".gemsembler"
        cache_path.mkdir(exist_ok=True, parents=True)
        cache_path /= func.__name__ + ".json"

        if cache_path.exists():
            with open(cache_path, "r") as fh:
                data = json.load(fh)
        else:
            data = func(*args, **kwargs)
            with open(cache_path, "w") as fh:
                json.dump(data, fh)
        return data

    return wrapper_decorator


def download_db(url, cache_name=None, **kwargs):
    """
    Function to download the data needed for conversion. Caches the data in
    ~/.gemsembler.
    """
    # Either get filename for cache file from input arg or url
    cache_name = cache_name or url.rsplit("/", 1)[-1]

    # Create the directory holding
    cache_path = Path("~").expanduser() / ".gemsembler"
    cache_path.mkdir(exist_ok=True, parents=True)
    cache_path /= cache_name

    # If the cache file exists open it, otherwise download the data
    if cache_path.exists():
        data = pd.read_csv(cache_path, sep="\t", **kwargs)
    else:
        data = pd.read_csv(url, sep="\t", **kwargs)
        data.to_csv(cache_path, sep="\t", index=False)
    return data


def process_bigg(data, metabolites=False):
    """
    Function to process the BiGG database and get the mapping between old and
    new BiGG ids.
    """
    return (
        data.dropna(subset="old_bigg_ids")
        .copy()
        .pipe(separate, "old_bigg_ids", sep="; ")
        .explode(column="old_bigg_ids")
        .get(["old_bigg_ids", "bigg_id"])
        .pipe(lambda x: replace(x, "bigg_id", "_[a-z]+$", "") if metabolites else x)
        .drop_duplicates()
        .groupby("old_bigg_ids", group_keys=False)
        .apply(lambda x: x["bigg_id"].tolist())
        .reset_index(name="bigg_ids")
        .rename(columns={"old_bigg_ids": "old_ids"})
        .set_index("old_ids")
        .get("bigg_ids")
        .to_dict()
    )


@cache_file
def get_bigg_m():
    df_bigg_m = download_db(
        "http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt",
        "bigg_models_metabolites.txt.gz",
    )
    return df_bigg_m.pipe(process_bigg, metabolites=True)


@cache_file
def get_bigg_r():
    df_bigg_r = download_db(
        "http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt",
        "bigg_models_reactions.txt.gz",
    )
    return df_bigg_r.pipe(process_bigg)


def process_modelseed(data):
    return (
        data[["id", "aliases"]]
        .set_index("id")["aliases"]
        .str.extract(r"BiGG: (.*?)\|")[0]
        .str.split("; ")
        .reset_index(name="bigg_ids")
        .rename(columns={"id": "seed_ids"})
        .set_index("seed_ids")
        .get("bigg_ids")
        .dropna()
        .to_dict()
    )


@cache_file
def get_seed_orig_m():
    df_modelseed_m = download_db(
        "https://github.com/ModelSEED/ModelSEEDDatabase/raw/master/Biochemistry/compounds.tsv",
        "compounds.tsv.gz",
    )
    return df_modelseed_m.pipe(process_modelseed)


@cache_file
def get_seed_orig_r():
    df_modelseed_r = download_db(
        "https://github.com/ModelSEED/ModelSEEDDatabase/raw/master/Biochemistry/reactions.tsv",
        "reactions.tsv.gz",
    )
    return df_modelseed_r.pipe(process_modelseed)


def process_metanetx(data, db_name, repl_regex):
    return (
        data.query("source.str.startswith(@db_name) or source.str.startswith('bigg')")
        .copy()
        .reset_index(drop=True)
        .pipe(colapply, "source", lambda x: x.str.replace(repl_regex, ":", regex=True))
        .pipe(separate, "source", sep=":", into=["db", "db_id"], expand=True)
        .drop_duplicates(subset=["db", "db_id"])
        .pivot_table(index="ID", columns="db", values="db_id", aggfunc=",".join)
        .dropna()
        .pipe(colapply, db_name, lambda x: x.str.split(","))
        .explode(db_name)
        .rename(columns={"bigg": "bigg_ids", db_name: f"{db_name}_ids"})
        .reset_index(drop=True)
        .set_index(f"{db_name}_ids")
        .get("bigg_ids")
        .str.split(",")
        .to_dict()
    )


@cache_file
def get_seed_addit_m():
    df_metanetx_m = download_db(
        "https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_xref.tsv",
        "chem_xref.tsv.gz",
        comment="#",
        names=["source", "ID", "description"],
    )
    return df_metanetx_m.pipe(
        process_metanetx, "seed", "(\.compound|M|\.metabolite):(M_)?"
    )


@cache_file
def get_seed_addit_r():
    df_metanetx_r = download_db(
        "https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_xref.tsv",
        "reac_xref.tsv.gz",
        comment="#",
        names=["source", "ID", "description"],
    )
    return df_metanetx_r.pipe(process_metanetx, "seed", "(\.reaction|R|):(R_)?")


@cache_file
def get_kegg_bigg_m():
    df_metanetx_m = download_db(
        "https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_xref.tsv",
        "chem_xref.tsv.gz",
        comment="#",
        names=["source", "ID", "description"],
    )
    return df_metanetx_m.pipe(
        process_metanetx,
        "kegg",
        "(\.compound|\.drug|\.metabolite|\.glycan|[CDGM]):(M_)?",
    )


@cache_file
def get_kegg_bigg_r():
    df_metanetx_r = download_db(
        "https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_xref.tsv",
        "reac_xref.tsv.gz",
        comment="#",
        names=["source", "ID", "description"],
    )
    return df_metanetx_r.pipe(process_metanetx, "kegg", "(\.reaction|R|):(R_)?")


def get_BiGG_lists(metabolites: bool):
    if metabolites:
        bigg_data = download_db(
            "http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt",
            "bigg_models_metabolites.txt.gz",
        )
        return set(bigg_data.get("universal_bigg_id"))
    else:
        bigg_data = download_db(
            "http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt",
            "bigg_models_reactions.txt.gz",
        )
        return set(bigg_data.get("bigg_id"))


def remove_duplicates(data, leave_from_mixed_directions=True):
    duplicated = (
        data.copy()
        .rename(
            columns={"1metabolites": "2metabolites", "2metabolites": "1metabolites"}
        )
        .pipe(lambda x: pd.concat([x, data], ignore_index=True))
        .sort_values(
            ["1metabolites", "2metabolites", "models_number", "reaction"],
            ascending=[True, True, False, True],
        )
        .pipe(lambda x: x[x.duplicated(["1metabolites", "2metabolites"], keep=False)])
    )

    reaction_list = duplicated["reaction"].unique().tolist()
    reaction_list_wo_dupl_metabolites = (
        duplicated.drop_duplicates(["1metabolites", "2metabolites"], keep="first")[
            "reaction"
        ]
        .unique()
        .tolist()
    )
    if leave_from_mixed_directions:
        reaction_list = set(reaction_list) - set(reaction_list_wo_dupl_metabolites)
    return data[~data["reaction"].isin(reaction_list)]


def get_bigg_network(path_to_dbs=None, leave_from_mixed_directions=True):
    """
    Getting dictionary for BiGG topology: unique reaction equation as key and
    id as value.  Maybe writing down a table with reaction ids unique reaction
    equations with metabolites sorted for the whole BiGG database. If reaction
    equations are duplicated one used in most amount of models selected.
    """

    bigg_database_r = download_db(
        "http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt",
        "bigg_models_reactions.txt.gz",
    )

    reaction_regex = re.compile(r"(\d+\.\d*(e-)?\d*|\d+e-\d*)|\+|<->")

    r_connections = (
        pd.DataFrame(
            {
                "reaction": bigg_database_r["bigg_id"],
                "models_number": bigg_database_r["model_list"].str.split().apply(len),
            }
        )
        # Add two columns for set of metabolites on both sides of the equation
        .pipe(
            set_value,
            ["1metabolites", "2metabolites"],
            (
                bigg_database_r["reaction_string"]
                .str.replace(r"(\d+\.\d*(e-)?\d*|\d+e-\d*)|\+", "", regex=True)
                .apply(lambda x: f" {x} ")
                .str.split("<->", n=1, expand=True)
            ),
        )
        # Sort the lists of metabolites in both columns
        .pipe(apply, "1metabolites", lambda x: " ".join(sorted(x.split())))
        .pipe(apply, "2metabolites", lambda x: " ".join(sorted(x.split())))
        # Sort the dataframe by the metabolites
        .sort_values(
            ["1metabolites", "2metabolites", "models_number"],
            ascending=[True, True, False],
        )
        # Remove duplicates based on metabolites
        .drop_duplicates(["1metabolites", "2metabolites"], keep="first")
        # Remove duplicates based on reactions
        .pipe(remove_duplicates, leave_from_mixed_directions)
        # Add an equation column
        .pipe(
            apply,
            ["1metabolites", "2metabolites"],
            lambda x: "<->".join(sorted(x)),
            "equation",
            axis=1,
        )
        .reset_index(drop=True)
    )

    return (
        r_connections.groupby("equation")
        .apply(lambda x: x["reaction"].values[0])
        .to_dict()
    )
