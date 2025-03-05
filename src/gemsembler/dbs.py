import json
from functools import wraps

import pandas as pd
from platformdirs import user_data_path


# helper functions for pandas dataframes
def separate(data, col, into=None, sep=" ", **kwargs):
    into = into or col
    new_data = data.copy()
    new_data[into] = data[col].str.split(sep, **kwargs)
    return new_data


def cache_file(func):
    """
    Decorator for caching mapping dictionaries in ~/.local/share/gemsembler
    """

    @wraps(func)
    def wrapper_decorator(*args, **kwargs):
        # Create the directory holding
        cache_path = user_data_path("gemsembler", ensure_exists=True)
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
    ~/.local/share/gemsembler.
    """
    # Either get filename for cache file from input arg or url
    cache_name = cache_name or url.rsplit("/", 1)[-1]

    # Define the path to where to cache the data
    cache_path = user_data_path("gemsembler", ensure_exists=True) / cache_name

    # If the cache file exists open it, otherwise download the data
    if cache_path.exists():
        data = pd.read_csv(cache_path, sep="\t", **kwargs)
    else:
        print(f"Downloading {cache_name} from {url}")
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
        .assign(
            bigg_id=lambda x: (
                x.bigg_id.str.replace("_[a-z]+$", "", regex=True)
                if metabolites
                else x.bigg_id
            )
        )
        .drop_duplicates()
        .groupby("old_bigg_ids", group_keys=False)["bigg_id"]
        .apply(list)
        .to_dict()
    )


@cache_file
def get_old_bigg_m():
    df_bigg_m = download_db(
        "http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt",
        "bigg_models_metabolites.txt.gz",
    )
    return df_bigg_m.pipe(process_bigg, metabolites=True)


@cache_file
def get_old_bigg_r():
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
        (
            "https://github.com/ModelSEED/ModelSEEDDatabase/"
            "raw/master/Biochemistry/compounds.tsv"
        ),
        "compounds.tsv.gz",
    )
    return df_modelseed_m.pipe(process_modelseed)


@cache_file
def get_seed_orig_r():
    df_modelseed_r = download_db(
        (
            "https://github.com/ModelSEED/ModelSEEDDatabase/"
            "raw/master/Biochemistry/reactions.tsv"
        ),
        "reactions.tsv.gz",
    )
    return df_modelseed_r.pipe(process_modelseed)


def process_with_metanetx(data, db_name, repl_regex):
    """Get a mapping between `db_name` and bigg using metanetx db."""
    return (
        data.query("source.str.startswith(@db_name) or source.str.startswith('bigg')")
        .copy()
        .reset_index(drop=True)
        .assign(source=lambda x: x.source.str.replace(repl_regex, ":", regex=True))
        .pipe(separate, "source", sep=":", into=["db", "db_id"], expand=True)
        .drop_duplicates(subset=["db", "db_id"])
        .pivot_table(index="ID", columns="db", values="db_id", aggfunc=",".join)
        .dropna()
        .assign(other_db=lambda x: x[db_name].str.split(","))
        .drop(columns=db_name)
        .explode("other_db")
        .set_index("other_db")
        .get("bigg")
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
        process_with_metanetx, "seed", r"(\.compound|M|\.metabolite):(M_)?"
    )


@cache_file
def get_seed_addit_r():
    df_metanetx_r = download_db(
        "https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_xref.tsv",
        "reac_xref.tsv.gz",
        comment="#",
        names=["source", "ID", "description"],
    )
    return df_metanetx_r.pipe(process_with_metanetx, "seed", r"(\.reaction|R|):(R_)?")


def process_metanetx(data, repl_regex):
    """Get a mapping between metanetx and bigg"""
    return (
        data.query("source.str.startswith('bigg') and ID != 'EMPTY'")
        .assign(source=lambda x: x.source.str.replace(repl_regex, "", regex=True))
        .rename(columns={"source": "bigg", "ID": "mnx"})
        .drop_duplicates(subset="bigg")
        .drop(columns="description")
        .groupby("mnx")["bigg"]
        .apply(list)
        .to_dict()
    )


@cache_file
def get_mnx_m():
    df_metanetx_m = download_db(
        "https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_xref.tsv",
        "chem_xref.tsv.gz",
        comment="#",
        names=["source", "ID", "description"],
    )
    return df_metanetx_m.pipe(
        process_metanetx, r"bigg(\.compound|M|\.metabolite):(M_)?"
    )


@cache_file
def get_mnx_r():
    df_metanetx_r = download_db(
        "https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_xref.tsv",
        "reac_xref.tsv.gz",
        comment="#",
        names=["source", "ID", "description"],
    )
    return df_metanetx_r.pipe(process_metanetx, r"bigg(\.reaction|R|):(R_)?")


@cache_file
def get_kegg_m():
    df_metanetx_m = download_db(
        "https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_xref.tsv",
        "chem_xref.tsv.gz",
        comment="#",
        names=["source", "ID", "description"],
    )
    return df_metanetx_m.pipe(
        process_with_metanetx,
        "kegg",
        r"(\.compound|\.drug|\.metabolite|\.glycan|[CDGM]):(M_)?",
    )


@cache_file
def get_kegg_r():
    df_metanetx_r = download_db(
        "https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_xref.tsv",
        "reac_xref.tsv.gz",
        comment="#",
        names=["source", "ID", "description"],
    )
    return df_metanetx_r.pipe(process_with_metanetx, "kegg", r"(\.reaction|R|):(R_)?")


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

    r_connections = (
        bigg_database_r["reaction_string"]
        .str.replace(r"(\d+\.\d*(e-)?\d*|\d+e-\d*)|\+", "", regex=True)
        .apply(lambda x: f" {x} ")
        .str.split("<->", n=1, expand=True)
        .rename(columns={0: "1metabolites", 1: "2metabolites"})
        # Sort the lists of metabolites in both columns
        .map(lambda x: " ".join(sorted(x.split())))
        .assign(
            reaction=bigg_database_r["bigg_id"],
            models_number=bigg_database_r["model_list"].str.split().apply(len),
        )
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
        .assign(
            equation=lambda x: (
                x[["1metabolites", "2metabolites"]].apply(
                    lambda row: "<->".join(sorted(row)), axis=1
                )
            )
        )
        .set_index("equation")
        .get("reaction")
        .to_dict()
    )

    return r_connections
